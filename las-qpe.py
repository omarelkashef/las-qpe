import numpy as np
from argparse import ArgumentParser
# PySCF imports
from pyscf import gto, scf, lib, mcscf, ao2mo
from pyscf.tools import fcidump
# mrh imports
from mrh.my_pyscf.mcscf.lasscf_o0 import LASSCF

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
from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper
from qiskit.providers.aer import StatevectorSimulator
from qiskit import Aer
from qiskit import QuantumCircuit
from qiskit.utils import QuantumInstance
from qiskit_nature.algorithms import VQEUCCFactory, GroundStateEigensolver
from qiskit_nature.circuit.library import HartreeFock
from qiskit.algorithms import NumPyEigensolver, PhaseEstimation, PhaseEstimationScale
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
    verbose=lib.logger.DEBUG)

# Do RHF
mf = scf.RHF(mol).run()
print("HF energy: ", mf.e_tot)

# Create LASSCF object
las = LASSCF(mf, (2,2),(2,2), spin_sub=(1,1))

# Localize the chosen fragment active spaces
frag_atom_list = ((0,1),(2,3))
loc_mo_coeff = las.localize_init_guess(frag_atom_list, mf.mo_coeff)

ncore = las.ncore
ncas = las.ncas
ncas_sub = las.ncas_sub

print ("Ncore: ",las.ncore, "Ncas: ",las.ncas, "Ncas_sub", las.ncas_sub)

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
D = mf.make_rdm1(mo_coeff=loc_mo_coeff)

hcore_ao = mol.intor_symmetric('int1e_kin') + mol.intor_symmetric('int1e_nuc')
hcore_mo = np.einsum('pi,pq,qj->ij', loc_mo_coeff, hcore_ao, loc_mo_coeff)

eri_4fold = ao2mo.kernel(mol.intor('int2e'), mo_coeffs=loc_mo_coeff)
eri = ao2mo.restore(1, eri_4fold,mol.nao_nr())

# Storing each fragment's h1 and h2 as a list
h1_frag = []
h2_frag = []

# Then construct h1' for each fragment
for idx in idx_list[1:]:
    inactive_mask = np.arange(mol.nao_nr())
    inactive_mask = np.delete(inactive_mask, idx)
    eri_mix = eri[:,:,:,inactive_mask]
    eri_mix = eri_mix[:,:,inactive_mask,:]
    h1p = hcore_mo[idx,idx]
    if h1p.size == 0:
        h1p = np.einsum('jkii->jk', eri_mix[idx,idx,:,:])
        for i,idx2 in enumerate(idx_list):
            if i > 0 and idx2 != idx:
                h1p = np.einsum('ijkl,kl->ij', eri[idx,idx,idx2,idx2],D[idx2,idx2])
    else:
        h1p += np.einsum('jkii->jk', eri_mix[idx,idx,:,:])
        for i,idx2 in enumerate(idx_list):
            if i > 0 and idx2 != idx:
                h1p += np.einsum('ijkl,kl->ij', eri[idx,idx,idx2,idx2],D[idx2,idx2])

    print("h1p is same as hcore: {}".format(np.allclose(h1p, hcore_mo)))
    # Finally, construct total H_frag
    h1_frag.append(h1p)
    h2_frag.append(eri[idx,idx,idx,idx])

print(h1_frag)
print(h2_frag)

for frag in range(len(ncas_sub)):
    # For QPE, need second_q_ops
    # Hacking together an ElectronicStructureDriverResult to create second_q_ops
    # Lines below stolen from qiskit's FCIDump driver and modified
    num_alpha = int(mol.nelec[0])
    num_beta = int(mol.nelec[1])

    particle_number = ParticleNumber(
        num_spin_orbitals=ncas_sub[frag]*2,
        num_particles=(num_alpha, num_beta),
    )

    # Assuming an RHF reference for now, so h1_b, h2_ab, h2_bb are None
    electronic_energy = ElectronicEnergy(
        [
            OneBodyElectronicIntegrals(ElectronicBasis.MO, (h1_frag[frag], None)),
            TwoBodyElectronicIntegrals(ElectronicBasis.MO, (h2_frag[frag], None, None, None)),
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
    qubit_converter = QubitConverter(mapper = ParityMapper(), two_qubit_reduction=True)
    # This just outputs a qubit op corresponding to a 2nd quantized op
    qubit_ops = [qubit_converter.convert(op) for op in second_q_ops]
    hamiltonian = qubit_ops[0]
    print(hamiltonian)

    # Set the backend
    quantum_instance = QuantumInstance(backend = Aer.get_backend('aer_simulator'), shots=args.shots)

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

    pe_scale = PhaseEstimationScale.from_pauli_sum(hamiltonian)

    # QK: scale so that phase does not wrap.
    scaled_hamiltonian = -pe_scale.scale * hamiltonian

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
    init_state = HartreeFock(8, (2,2), qubit_converter)

    # Estimate takes in a SummedPauli or a PauliOp and outputs a scaled estimate of the eigenvalue
    res = qpe_solver.estimate(unitary=unitary, state_preparation=init_state)

    print(res)
    scaled_phases = pe_scale.scale_phases(res.filter_phases(cutoff=0.0,as_float=True))
    print(scaled_phases)

las.kernel(loc_mo_coeff)
print("LASSCF energy: ", las.e_tot)
