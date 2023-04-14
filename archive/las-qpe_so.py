import numpy as np
import logging
from argparse import ArgumentParser
# PySCF imports
from pyscf import gto, scf, lib, mcscf, ao2mo
from pyscf.tools import fcidump
# mrh imports
from mrh.my_pyscf.mcscf.lasscf_o0 import LASSCF
from mrh.my_pyscf.mcscf.lasci import h1e_for_cas
from c4h6_struct import structure

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
from qiskit.quantum_info import DensityMatrix, partial_trace
from qiskit.quantum_info.operators.channel import SuperOp
from qiskit.utils import QuantumInstance
from qiskit_nature.algorithms import VQEUCCFactory, GroundStateEigensolver
from qiskit_nature.circuit.library import HartreeFock, UCCSD
from qiskit.algorithms import NumPyEigensolver, PhaseEstimation, PhaseEstimationScale, VQE 
from qiskit.algorithms.optimizers import L_BFGS_B
from qiskit.algorithms.phase_estimators import PhaseEstimationResult
from qiskit.opflow import PauliTrotterEvolution,SummedOp,PauliOp,MatrixOp,PauliSumOp,StateFn

parser = ArgumentParser(description='Do QPE with a LAS reference')
parser.add_argument('--an', type=int, default=1, help='number of ancilla qubits')
parser.add_argument('--shots', type=int, default=1024, help='number of shots for the simulator')
args = parser.parse_args()

# Define molecule: (H2)_2
xyz = '''H 0.0 0.0 0.0
         H 1.0 0.0 0.0
         H 0.2 3.9 0.1
        H 1.159166 4.1 -0.1'''
mol = gto.M (atom = xyz, basis = 'sto-3g', output='h4_sto3g.log',
    symmetry=False, verbose=lib.logger.DEBUG)

# Define molecule: C_4H_6
#norb = 8
#nelec = 8
#norb_f = (4,4)
#nelec_f = ((2,2),(2,2))
#mol = structure (0.0, 0.0, output='c4h6_631g_{}_{}.log'.format(args.an, args.shots), verbose=lib.logger.DEBUG)

# Do RHF
mf = scf.RHF(mol).run()
print("HF energy: ", mf.e_tot)

# Create LASSCF object
#las = LASSCF(mf, (2,),(2,), spin_sub=(1,))
#las = LASSCF(mf, (1,1),(1,1), spin_sub=(2,2))
las = LASSCF(mf, (2,2),(2,2), spin_sub=(1,1))
#las = LASSCF(mf, (4,),(4,), spin_sub=(1,))

# Localize the chosen fragment active spaces
#frag_atom_list = ((0,1),)
#frag_atom_list = ((0,),(1,))
frag_atom_list = ((0,1),(2,3))
#frag_atom_list = ((0,1,2,3),)
loc_mo_coeff = las.localize_init_guess(frag_atom_list, mf.mo_coeff)
las.kernel(loc_mo_coeff)
loc_mo_coeff = las.mo_coeff
print("LASSCF energy: ", las.e_tot)
#loc_mo_coeff = mf.mo_coeff

ncore = las.ncore
ncas = las.ncas
ncas_sub = las.ncas_sub
nelec_cas = las.nelecas

print ("Ncore: ", ncore, "Ncas: ", ncas, "Ncas_sub: ", ncas_sub, "Nelec_cas: ", nelec_cas)

# Situation so far: we have mf.mo_coeffs containing [:,ncore:nsub1:nsub2:next]
# Creating a list of slices for core, subspace1, subspace2, etc
idx_list = [slice(0,ncore)]
prev_sub_size = 0
for i, sub in enumerate(ncas_sub):
    idx_list.append(slice(ncore+prev_sub_size, ncore+prev_sub_size + sub))
    prev_sub_size += sub

# To prepare: H_eff = \sum_K H_frag(K)
# H_frag(K) = h1'_k1^{k2} a^+_{k1} a_{k2} + 1/4 h2_{k2 k4}^{k1 k3} a^+_{k1} a^+_{k3} a_{k4} a_{k2}
# with h1' = h1_{k1}^{k2} + \sum_i h2_{k2 i}^{k1 i} + \sum{L \neq K} h2_{k2 l2}^{k1 l1} D_{l2}^{l1}

# First, construct D and ints
# Option 1: AO-basis HF 1-RDM to localized MO basis
'''
D = mf.make_rdm1(mo_coeff=mf.mo_coeff)
D_mo = np.einsum('pi,pq,qj->ij', loc_mo_coeff, D, loc_mo_coeff)
'''
# Option 2: Converged LAS 1-RDM
D_mo = las.make_rdm1(mo_coeff=las.mo_coeff)
# Convert to spin orbitals
D_so = np.block([[D_mo, np.zeros_like(D_mo)],[np.zeros_like(D_mo), D_mo]])
#print("D:\n",D_so)

hcore_ao = mf.get_hcore(mol)
hcore_mo = np.einsum('pi,pq,qj->ij', loc_mo_coeff, hcore_ao, loc_mo_coeff)
# Convert to spin orbitals
hcore_so = np.block([[hcore_mo, np.zeros_like(hcore_mo)],[np.zeros_like(hcore_mo), hcore_mo]])

nso = mol.nao_nr() * 2
eri_4fold = ao2mo.kernel(mol.intor('int2e'), mo_coeffs=loc_mo_coeff)
eri = ao2mo.restore(1, eri_4fold,mol.nao_nr())
#print(eri[0])
# Convert to spin orbitals
eri_so = np.zeros((nso, nso, nso, nso))
one_indices = (
            (0, 0, 0, 0),  # alpha-alpha-spin
            (0, 1, 1, 0),  # beta-alpha-spin
            (1, 1, 1, 1),  # beta-beta-spin
            (1, 0, 0, 1),  # alpha-beta-spin
        )
for one_idx in one_indices: 
    ao_mat = np.einsum('ijkl->ljik', eri)  # physicist notation
    kron = np.zeros((2,2,2,2))
    kron[one_idx] = 1
    eri_so -= 0.5 * np.kron(kron, ao_mat)
#print(eri_so[0])

# Storing each fragment's h1 and h2 as a list
h1_frag = []
h2_frag = []

# Then construct h1' for each fragment
# and for alpha-alpha, beta-beta blocks
for idx in idx_list[1:]:
    print(idx_list[0])
    inactive_mask = idx_list[0]
    eri_mix = eri_so[:,inactive_mask,:,inactive_mask].copy()
    print(eri_mix.shape)
    h1p = hcore_so[idx,idx].copy()
    if h1p.size == 0:
        h1p = np.einsum('jiki->jk', eri_mix[idx,:,idx,:])
        for i,idx2 in enumerate(idx_list):
            if i > 0 and idx2 != idx:
                h1p += np.einsum('ikjl,kl->ij', eri_so[idx,idx2,idx,idx2],D_so[idx2,idx2])
    else:
        h1p += np.einsum('jiki->jk', eri_mix[idx,:,idx,:])
        for i,idx2 in enumerate(idx_list):
            if i > 0 and idx2 != idx:
                h1p += np.einsum('ikjl,kl->ij', eri_so[idx,idx2,idx,idx2],D_so[idx2,idx2])

    # Finally, construct total H_frag
    h1_frag.append(h1p)
    h2_frag.append(eri[idx,idx,idx,idx])

# using the built-in LASCI function h1e_for_cas
h1_frag_2 = h1e_for_cas(las)

# Just using h1e_for_cas as my fragment h1
h1_frag = []
for f in range(len(ncas_sub)):
    h1_frag.append(h1_frag_2[f][0][0])
print("h1e: {}".format(h1_frag))

for f in range(len(ncas_sub)):
    print("H1_frag shape: ", h1_frag[f].shape)
    print("H2_frag shape: ", h2_frag[f].shape)

# Function below stolen from qiskit's Hamiltonian Phase Estimation class
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
dm_list = []

for frag in range(len(ncas_sub)):
    # WARNING: these have to be set manually for each fragment!
    num_alpha = int(nelec_cas[0] / 2)
    num_beta = int(nelec_cas[1] / 2)

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
    second_q_ops[0].set_truncation(0)
    #print(second_q_ops[0])

    # Choose fermion-to-qubit mapping
    qubit_converter = QubitConverter(mapper = ParityMapper(), two_qubit_reduction=True)
    # This just outputs a qubit op corresponding to a 2nd quantized op
    qubit_ops = [qubit_converter.convert(op) for op in second_q_ops]
    hamiltonian = qubit_ops[0]
    #print(hamiltonian)

    # Set the backend
    quantum_instance = QuantumInstance(backend = Aer.get_backend('qasm_simulator'), shots=args.shots, optimization_level=0)

    np_solver = NumPyEigensolver(k=1)
    ed_result = np_solver.compute_eigenvalues(hamiltonian)
    print("NumPy result: ", ed_result.eigenvalues)

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

    pe_scale = PhaseEstimationScale.from_pauli_sum(hamiltonian_no_id)

    # QK: scale so that phase does not wrap.
    scaled_hamiltonian = -pe_scale.scale * hamiltonian_no_id  

    # Default evolution: PauliTrotterEvolution
    evolution = PauliTrotterEvolution()

    # Create the unitary by evolving the Hamiltonian
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
    # For our (H2)_2 system 8 spin orbs, 2 alpha 2 beta electrons
    init_state = HartreeFock(ncas_sub[frag]*2, (num_alpha,num_beta), qubit_converter)

    circuit = qpe_solver.construct_circuit(unitary=unitary, state_preparation=init_state).decompose().decompose()

    op_dict = circuit.count_ops()
    total_ops = sum(op_dict.values())
    total_op_list.append(total_ops)
    print("Operations: {}".format(op_dict))
    print("Total operations: {}".format(total_ops))
    print("Nonlocal gates: {}".format(circuit.num_nonlocal_gates()))

    # Estimate takes in a SummedPauli or a PauliOp and outputs a scaled estimate of the eigenvalue
    res = qpe_solver.estimate(unitary=unitary, state_preparation=init_state)
    phases = res.__dict__['_phases']
    phases_list.append(phases)
    print(res)
    scaled_phases = pe_scale.scale_phases(res.filter_phases(cutoff=0.0, as_float=True), id_coefficient=0.0)
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
    max_count = 5
    count = 0
    while np.allclose(new_eig, most_likely_eig) is False:
        print("Reusing... [",count,"]")
        # Generating a new instance
        new_instance = QuantumInstance(backend = Aer.get_backend('qasm_simulator'), shots=1)

        # Using the new instance in a solver
        new_qpe_solver = PhaseEstimation(num_evaluation_qubits=args.an, quantum_instance=new_instance)

        # Creating the circuit in order to use the save_statevector instruction
        new_circuit = QuantumCircuit(ncas_sub[frag]*2+args.an)
        print(new_circuit.draw())
        new_circuit = new_qpe_solver.construct_circuit(unitary=unitary, state_preparation=init_state)
        # Reusing the already-prepared unitary and initial state
        ## To obtain a statevector after measurement, I must use the class function
        ## to add the measurements into the circuit before appending the save_statevector
        ## instruction. This is ugly, as I'm accessing class functions not meant to be
        ## directly accessed.
        new_qpe_solver._add_measurement_if_required(new_circuit)
        new_circuit.save_statevector(label='final')
        print(new_circuit.decompose().draw())
        circuit_result = new_qpe_solver._quantum_instance.execute(new_circuit)
        phases = new_qpe_solver._compute_phases(ncas_sub[frag]*2, circuit_result)
        gs_result = PhaseEstimationResult(args.an, circuit_result=circuit_result, phases=phases)
        pe_scale = PhaseEstimationScale.from_pauli_sum(hamiltonian_no_id)
        scaled_phases = pe_scale.scale_phases(gs_result.filter_phases(cutoff=0.0, as_float=True), id_coefficient=0.0)
        (new_eig, v), = scaled_phases.items()
        print("New eig: ",new_eig)

        # Saving a density matrix as an alternative to saving a statevector
        # Can't load density matrices on 2 different registers on a single circuit 
        # Also this density matrix is created without the measurement gates
        #DM = DensityMatrix.from_instruction(new_qpe_solver.construct_circuit(unitary=unitary, state_preparation=init_state))
        #print("Density matrix shape: ", DM.data.shape)
        #print("Density matrix purity: ", DM.purity())
        #eigvals, eigvecs = np.linalg.eig(DM.data)
        #print(eigvals)
        #print("Statevector:",eigvecs)
        #PT = partial_trace(DM, list(range(args.an)))
        #print("Partial trace purity: ", PT.purity())
        #psi = gs_result.circuit_result.data(0)['final']

        # Checking that partial_trace does what it's supposed to do
        #state = DM
        #qargs = list(range(args.an))
        ### Taken from qiskit.quantum_info.partial_trace
        ## QK:Compute traced shape
        #traced_shape = state._op_shape.remove(qargs=qargs)
        #print("Traced_shape",traced_shape)
        #print("Traced_shape qargs",traced_shape._num_qargs_l)
        ## QK:Density matrix case
        ## QK:Trace first subsystem to avoid coping whole density matrix
        #dims = state.dims(qargs)
        #print("dims",dims)
        #tr_op = SuperOp(np.eye(dims[0]).reshape(1, dims[0] ** 2), input_dims=[dims[0]], output_dims=[1])
        #ret = state.evolve(tr_op, [qargs[0]])
        ## QK:Trace over remaining subsystems
        #for qarg, dim in zip(qargs[1:], dims[1:]):
        #    tr_op = SuperOp(np.eye(dim).reshape(1, dim ** 2), input_dims=[dim], output_dims=[1])
        #    ret = ret.evolve(tr_op, [qarg])
        ## QK:Remove traced over subsystems which are listed as dimension 1
        #ret._op_shape = traced_shape

        count = count + 1
        if count > max_count:
            print("Max iterations exceeded.")
            break

    # Save only the statevector corresponding to the system
    final_wfn = gs_result.circuit_result.data(0)['final']
    (an_state, v), = gs_result.__dict__['_phases'].items()
    print("Ancilla state: ",int(an_state,2))
    #print(final_wfn)
    final_wfn = final_wfn[int(an_state, 2)::2**args.an]
    #print(final_wfn)
    state_list.append(final_wfn)
    #print("Partial trace shape: ",PT.data.shape)
    #dm_list.append(PT)

print("Phases: ",phases_list)
print("en_list: ",en_list)
print("total_op_list: ",total_op_list)

# Setting up the Hamiltonian for the 2-fragment system
driver = PySCFDriver(atom=xyz, charge=0, spin=0, method=MethodType.RHF)
driver_result = driver.run()

second_q_ops = driver_result.second_q_ops()
# This just outputs a qubit op corresponding to a 2nd quantized op
qubit_ops = [qubit_converter.convert(op) for op in second_q_ops]
hamiltonian = qubit_ops[0]

# Create a quantum register with system qubits
qr1 = QuantumRegister(np.sum(ncas_sub)*2, 'q1')

# Create the state by  and normalizing
total_state = np.kron(state_list[0],state_list[1])
total_state = total_state/np.linalg.norm(total_state)

new_circuit = QuantumCircuit(qr1)
new_circuit.initialize(total_state, qubits=qr1)
#print(new_circuit.draw())

# Need to set up a new qubit_converter and quantum instance
qubit_converter = QubitConverter(mapper = ParityMapper(), two_qubit_reduction=True)
new_instance = QuantumInstance(backend = Aer.get_backend('aer_simulator'), shots=1024)

# Setting up the VQE
ansatz = UCCSD(qubit_converter=qubit_converter, num_particles=(2,2), num_spin_orbitals=8, initial_state=new_circuit)
optimizer = L_BFGS_B(maxfun=10000, iprint=101)

algorithm = VQE(ansatz=ansatz, optimizer=optimizer, quantum_instance=new_instance) 
result = algorithm.compute_minimum_eigenvalue(hamiltonian)
print(result)

np.save('results_{}_{}.npy'.format(args.an, args.shots), {'n_frag':len(ncas_sub), 'phases':phases_list, 'energies':en_list, 'operations':total_op_list})

