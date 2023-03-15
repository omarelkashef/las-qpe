#########################
# Script to run LASSCF, use the converged fragment Hamiltonians to set up
# fragment wavefunctions using QPE, then save the state 
#########################

import numpy as np
import logging
import time
from argparse import ArgumentParser
# PySCF imports
from pyscf import gto, scf, lib, mcscf, ao2mo
from pyscf.tools import fcidump
# mrh imports
from mrh.my_pyscf.mcscf.lasscf_o0 import LASSCF
from mrh.my_pyscf.mcscf.lasci import h1e_for_cas
#from c4h6_struct import structure
from get_geom import get_geom

# Qiskit imports
from qiskit_nature.properties.second_quantization.electronic import (
    ElectronicStructureDriverResult,
    ElectronicEnergy,
    ParticleNumber,
)
from qiskit_nature.properties.second_quantization.electronic.integrals import (
    OneBodyElectronicIntegrals,
    TwoBodyElectronicIntegrals,
)
from qiskit_nature.properties.second_quantization.electronic.bases import ElectronicBasis
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.drivers.second_quantization import PySCFDriver, MethodType
from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper
from qiskit.providers.aer import StatevectorSimulator, QasmSimulator
from qiskit import Aer, transpile
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
#from qiskit.quantum_info.states.densitymatrix import DensityMatrix
from qiskit.visualization import plot_state_city
from qiskit.quantum_info import DensityMatrix, partial_trace, Statevector
from qiskit.quantum_info.operators.channel import SuperOp
from qiskit.utils import QuantumInstance
from qiskit_nature.algorithms import VQEUCCFactory, GroundStateEigensolver
from qiskit_nature.circuit.library import HartreeFock, UCCSD, UCC
from qiskit.algorithms import NumPyEigensolver, PhaseEstimation, PhaseEstimationScale, VQE 
from qiskit.algorithms.optimizers import L_BFGS_B, COBYLA, BOBYQA
from qiskit.algorithms.phase_estimators import PhaseEstimationResult
from qiskit.opflow import PauliTrotterEvolution,SummedOp,PauliOp,MatrixOp,PauliSumOp,StateFn

parser = ArgumentParser(description='Do LAS-QPE, specifying num of ancillas and shots')
parser.add_argument('--an', type=int, default=1, help='number of ancilla qubits')
parser.add_argument('--dist', type=float, default=1.0, help='distance of H2s from one another')
parser.add_argument('--shots', type=int, default=1024, help='number of shots for the simulator')
args = parser.parse_args()

# Define molecule: (H2)_3
xyz = get_geom('scan', dist=args.dist)
mol = gto.M (atom = xyz, basis = 'sto-3g', output='h4_sto3g_{}.log'.format(args.dist),
    symmetry=False, verbose=lib.logger.DEBUG)

# Do RHF
mf = scf.RHF(mol).run()
print("HF energy: ", mf.e_tot)

# Create LASSCF object
# Keywords: (wavefunction obj, num_orb in each subspace, (num_alpha in each subspace, num_beta in each subspace), spin multiplicity in each subspace)
#las = LASSCF(mf, (2,),(2,), spin_sub=(1,))
las = LASSCF(mf, (2,2,2),(2,2,2), spin_sub=(1,1,1))

# Localize the chosen fragment active spaces
frag_atom_list = ((0,1),(2,3),(4,5))
loc_mo_coeff = las.localize_init_guess(frag_atom_list, mf.mo_coeff)

# Run LASSCF
las.kernel(loc_mo_coeff)
loc_mo_coeff = las.mo_coeff
print("LASSCF energy: ", las.e_tot)

ncore = las.ncore
ncas = las.ncas
ncas_sub = las.ncas_sub
nelecas_sub = las.nelecas_sub

print ("Ncore: ", ncore, "Ncas: ", ncas, "Ncas_sub: ", ncas_sub, "Nelecas_sub: ", nelecas_sub)

# Situation so far: we have loc_mo_coeff containing [:,ncore:nsub1:nsub2:next]
# Creating a list of slices for core, subspace1, subspace2, etc
idx_list = [slice(0,ncore)]
prev_sub_size = 0
for i, sub in enumerate(ncas_sub):
    idx_list.append(slice(ncore+prev_sub_size, ncore+prev_sub_size + sub))
    prev_sub_size += sub

# using the built-in LASCI functions h1e_for_cas, get_h2eff
h1_las = las.h1e_for_cas()
eri_las = las.get_h2eff(loc_mo_coeff)

# Storing each fragment's h1 and h2 as a list
h1_frag = []
h2_frag = []

# Then construct h1' for each fragment
# and for alpha-alpha, beta-beta blocks
for idx in range(len(ncas_sub)):
    h2_frag.append(las.get_h2eff_slice(eri_cas, idx))

# Just using h1e_for_cas as my fragment h1
h1_frag = []
for f in range(len(ncas_sub)):
    h1_frag.append(h1_las[f][0][0])

# Checking that the fragment Hamiltonian shapes are correct
for f in range(len(ncas_sub)):
    print("H1_frag shape: ", h1_frag[f].shape)
    print("H2_frag shape: ", h2_frag[f].shape)

# Function below stolen from qiskit's Hamiltonian Phase Estimation class
# To make the QPE slightly less expensive
def _remove_identity(pauli_sum):
    """Remove any identity operators from `pauli_sum`. Return
    the sum of the coefficients of the identities and the new operator.
    """
    idcoeff = 0.0
    ops = []
    for op in pauli_sum:
        p = op.primitive
        if p.x.any() or p.z.any():
            ops.append(op)
        else:
            idcoeff += op.coeff

    return idcoeff, SummedOp(ops)

phases_list = []
en_list = []
total_op_list = []
state_list = []
result_list = []

frag_t0 = time.time()

for frag in range(len(ncas_sub)):
    # Get alpha and beta electrons from LAS
    num_alpha = nelecas_sub[frag][0]
    num_beta = nelecas_sub[frag][1]

    # For QPE, need second_q_ops
    # Hacking together an ElectronicStructureDriverResult to create second_q_ops
    # Lines below stolen from qiskit's FCIDump driver and modified
    particle_number = ParticleNumber(
        num_spin_orbitals=ncas_sub[frag]*2,
        num_particles=(num_alpha, num_beta),
    )

    # Assuming an RHF reference for now, so h1_b, h2_ab, h2_bb are created using 
    # the corresponding spots from h1_frag and just the aa term from h2_frag
    print("Nuclear repulsion: ", las.energy_nuc())
    electronic_energy = ElectronicEnergy(
        [
            # Using MO basis here for simplified conversion
            OneBodyElectronicIntegrals(ElectronicBasis.MO, (h1_frag[frag], None)),
            TwoBodyElectronicIntegrals(ElectronicBasis.MO, (h2_frag[frag], h2_frag[frag], h2_frag[frag], None)),
        ],
        nuclear_repulsion_energy=las.energy_nuc(),
    )

    # QK NOTE: under Python 3.6, pylint appears to be unable to properly identify this case of
    # nested abstract classes (cf. https://github.com/Qiskit/qiskit-nature/runs/3245395353).
    # However, since the tests pass I am adding an exception for this specific case.
    # pylint: disable=abstract-class-instantiated
    driver_result = ElectronicStructureDriverResult()
    driver_result.add_property(electronic_energy)
    driver_result.add_property(particle_number)

    second_q_ops = driver_result.second_q_ops()

    # Choose fermion-to-qubit mapping
    qubit_converter = QubitConverter(mapper = JordanWignerMapper(), two_qubit_reduction=False)
    # This just outputs a qubit op corresponding to a 2nd quantized op
    qubit_ops = [qubit_converter.convert(op) for op in second_q_ops]
    hamiltonian = qubit_ops[0]

    # Set the backend
    quantum_instance = QuantumInstance(backend = Aer.get_backend('aer_simulator'), shots=args.shots)

    # Numpy solver to estimate error in QPE energy due to trotterization
    np_solver = NumPyEigensolver(k=1)
    ed_result = np_solver.compute_eigenvalues(hamiltonian)
    print("NumPy result: ", ed_result.eigenvalues)
    numpy_wfn = ed_result.eigenstates

    # Can choose a regular solver from qiskit.algorithms
    qpe_solver = PhaseEstimation(num_evaluation_qubits=args.an, quantum_instance=quantum_instance)

    # Create a unitary out of the Hamiltonian ops
    # Lines below stolen from qiskit.algorithms.HPE
    if isinstance(hamiltonian, PauliSumOp):
        hamiltonian = hamiltonian.to_pauli_op()
    elif isinstance(hamiltonian, PauliOp):
        hamiltonian = SummedOp([hamiltonian])

    if isinstance(hamiltonian, SummedOp):
        id_coefficient, hamiltonian_no_id = _remove_identity(hamiltonian)
    else:
        raise TypeError("Hamiltonian must be PauliSumOp, PauliOp or SummedOp.")

    # Instantiate a PEScale object for conversion later
    pe_scale = PhaseEstimationScale.from_pauli_sum(hamiltonian_no_id)

    # QK: scale so that phase does not wrap.
    scaled_hamiltonian = -pe_scale.scale * hamiltonian_no_id  

    # Default evolution: PauliTrotterEvolution
    evolution = PauliTrotterEvolution(reps=7)

    # Create the unitary by evolving the Hamiltonian
    # Here is the source of Trotter error
    unitary = evolution.convert(scaled_hamiltonian.exp_i())

    if not isinstance(unitary, QuantumCircuit):
        unitary_circuit = unitary.to_circuit()
    else:
        unitary_circuit = unitary

    # QK: Decomposing twice allows some 1Q Hamiltonians to give correct results
    # QK: when using MatrixEvolution(), that otherwise would give incorrect results.
    # QK: It does not break any others that we tested.
    unitary = unitary_circuit.decompose().decompose()

    # Printing this is not a good idea because the circuit is very large
    #print(unitary)

    # Create an HF initial state and add it to the estimate function
    # For our H_2 system, 4 spin orbs, 1 alpha 1 beta electron
    init_state = HartreeFock(ncas_sub[frag]*2, (num_alpha,num_beta), qubit_converter)

    # Gate counts
    if int(args.dist) == 0.0:
        circuit = qpe_solver.construct_circuit(unitary=unitary, state_preparation=init_state).decompose()
        target_basis = ['rx', 'ry', 'rz', 'h', 'cx']
        circ_for_counts = transpile(circuit, basis_gates=target_basis, optimization_level=0)

        op_dict = circ_for_counts.count_ops()
        total_ops = sum(op_dict.values())
        total_op_list.append(total_ops)
        print("Operations: {}".format(op_dict))
        print("Total operations: {}".format(total_ops))
        print("Nonlocal gates: {}".format(circuit.num_nonlocal_gates()))

    # Estimate takes in a SummedPauli or a PauliOp and outputs a scaled estimate of the eigenvalue
    frag_result = qpe_solver.estimate(unitary=unitary, state_preparation=init_state)
    phases = frag_result.__dict__['_phases']
    phases_list.append(phases)
    print(frag_result)
    scaled_phases = pe_scale.scale_phases(frag_result.filter_phases(cutoff=0.0, as_float=True), id_coefficient=id_coefficient)
    scaled_phases = {v:k for k, v in scaled_phases.items()}
    print(scaled_phases)
    energy_dict = {k:scaled_phases[v] for k, v in phases.items()}
    en_list.append(energy_dict)
    most_likely_eig = scaled_phases[max(scaled_phases.keys())]
    most_likely_an = max(phases, key=phases.get)
    print("Most likely eigenvalue: ", most_likely_eig)
    print("Most likely ancilla sign: ", most_likely_an)

    # For a given fragment, rerun the QPE until you get the ground state

    new_eig = 1e-5
    max_count = 10
    count = 0
    while np.allclose(new_eig, most_likely_eig) is False:
        print("Reusing... [",count,"]")
        # Generating a new single-shot instance
        new_instance = QuantumInstance(backend = Aer.get_backend('aer_simulator'), shots=1)

        # Using the new instance in a solver
        new_qpe_solver = PhaseEstimation(num_evaluation_qubits=args.an, quantum_instance=new_instance)

        # Reusing the already-prepared unitary and initial state
        new_circuit = new_qpe_solver.construct_circuit(unitary=unitary, state_preparation=init_state)

        ## To obtain a statevector after measurement, I must use the class function
        ## to add the measurements into the circuit before appending the save_statevector
        ## instruction. This is ugly, as I'm accessing class functions not meant to be
        ## directly accessed, but necessary.
        new_qpe_solver._add_measurement_if_required(new_circuit)
        new_circuit.save_statevector(label='final')

        # Run the circuit with the save instruction
        circuit_result = new_qpe_solver._quantum_instance.execute(new_circuit)
        phases = new_qpe_solver._compute_phases(ncas_sub[frag]*2, circuit_result)
        gs_result = PhaseEstimationResult(args.an, circuit_result=circuit_result, phases=phases)
        pe_scale = PhaseEstimationScale.from_pauli_sum(hamiltonian_no_id)
        scaled_phases = pe_scale.scale_phases(gs_result.filter_phases(cutoff=0.0, as_float=True), id_coefficient=id_coefficient)
        (new_eig, v), = scaled_phases.items()
        print("New eig: ",new_eig)

        count = count + 1
        if count > max_count:
            print("Max iterations exceeded.")
            break

    # Save only the statevector corresponding to the system qubits
    final_wfn = gs_result.circuit_result.data(0)['final']
    (an_state, v), = gs_result.__dict__['_phases'].items()
    print("Ancilla state: ",an_state)
    print("Before reducing:",final_wfn)
    final_wfn = final_wfn._data[int(an_state[::-1],2)::2**args.an]
    print("After reducing: ",final_wfn)
    final_state = Statevector(final_wfn)
    overlap = numpy_wfn[0].primitive.inner(final_state)
    print("Overlap of numpy wfn and QPE statevector: ", overlap)
    state_list.append(final_wfn)
    result_list.append(gs_result)

frag_t1 = time.time()

print("Fragment QPE total time (s): ",frag_t1-frag_t0)
print("Phases: ",phases_list)
print("en_list: ",en_list)
print("total_op_list: ",total_op_list)

np.save("qpe_state_{}.npy".format(args.dist), state_list, allow_pickle=True)
