##########################
# Script to run LASSCF 
# and use a LASCI vector initialized on the qubits
# to run a VQE using a generalized UCCSD ansatz on the whole system
# It can also use an HF initial state
#########################

import numpy as np
import logging
import time
from argparse import ArgumentParser
from typing import Tuple, List
import itertools
# PySCF imports
from pyscf import gto, scf, lib, mcscf, ao2mo
from pyscf.tools import fcidump
# mrh imports
from mrh.my_pyscf.mcscf.lasscf_o0 import LASSCF
from mrh.my_pyscf.mcscf.lasci import h1e_for_cas
from mrh.exploratory.unitary_cc import uccsd_sym1, lasuccsd
#from c4h6_struct import structure
from get_geom import get_geom
from custom_UCC import custom_UCC
from get_hamiltonian import get_hamiltonian

# Qiskit imports
#from qiskit_nature.converters.second_quantization import QubitConverteri
#from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper
from qiskit_nature.second_q.mappers import ParityMapper , QubitConverter , JordanWignerMapper
#from qiskit.providers.aer import StatevectorSimulator, QasmSimulator
#from qiskit import Aer, transpile
from qiskit.compiler import transpile
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
#from qiskit.quantum_info import DensityMatrix, partial_trace, Statevector
#from qiskit.utils import QuantumInstance
#from qiskit_nature.algorithms import VQEUCCFactory, GroundStateEigensolver
#from qiskit_nature.circuit.library import HartreeFock, UCCSD
from qiskit_nature.second_q.circuit.library import UCCSD , HartreeFock
#from qiskit.algorithms import NumPyEigensolver, VQE 
from qiskit.algorithms.minimum_eigensolvers import VQE 
from qiskit.algorithms.optimizers import L_BFGS_B, COBYLA, BOBYQA
from qiskit.opflow import PauliTrotterEvolution,SummedOp,PauliOp,MatrixOp,PauliSumOp,StateFn
from qiskit.primitives import Estimator
from qiskit import Aer




parser = ArgumentParser(description='Do LAS-VQE, specifying num of ancillas and shots')
parser.add_argument('--dist', type=float, default=1.35296239, help='distance of H2s from one another or C=C bond distance scaling in butadiene')
parser.add_argument('--shots', type=int, default=1024, help='number of shots for the simulator')
args = parser.parse_args()

# Prints info at every VQE iteration
logging.basicConfig(level='DEBUG')

# Define molecule: (H2)_2
xyz = get_geom('scan', dist=args.dist)
#xyz = '''H 0.0 0.0 0.0
#             H 1.0 0.0 0.0
#             H 0.2 1.6 0.1
#             H 1.159166 1.3 -0.1'''
mol = gto.M (atom = xyz, basis = 'sto-3g', output='h4_sto3g_{}.log'.format(args.dist),
    symmetry=False, verbose=lib.logger.DEBUG)

# Do RHF
mf = scf.RHF(mol).run()
print("HF energy: ", mf.e_tot)

'''
# Define molecule: C_4H_6
norb = 8
nelec = 8
norb_f = (4,4)
nelec_f = ((2,2),(2,2))
mol = structure (0.0, 0.0, output='c4h6_631g_{}.log'.format(args.dist), verbose=lib.logger.DEBUG)

# Do RHF
mf = scf.RHF(mol).run()
print("HF energy: ", mf.e_tot)

#CASSCF (for comparison with Matt's code)
mc = mcscf.CASSCF (mf, norb, nelec).set (fcisolver = csf_solver (mol, smult=1))
mo_coeff = mc.sort_mo ([11,12,14,15,16,17,21,24])
mc.kernel (mo_coeff)
mo_coeff = mc.mo_coeff.copy ()

# Create LASSCF object
las = LASSCF (mf, norb_f, nelec_f, spin_sub=(1,1)).set (mo_coeff=mo_coeff)

# Localize the chosen fragment active spaces
loc_mo_coeff = las.localize_init_guess ([[0,1,2,3,4],[5,6,7,8,9]])
las.kernel(loc_mo_coeff)
loc_mo_coeff = las.mo_coeff
print("LASSCF energy: ", las.e_tot)
'''
# Create LASSCF object
# Keywords: (wavefunction obj, num_orb in each subspace, (nelec in each subspace)/((num_alpha, num_beta) in each subspace), spin multiplicity in each subspace)
#las = LASSCF(mf, (2,),(2,), spin_sub=(1,))
las = LASSCF(mf, (2,2),(2,2), spin_sub=(1,1))
#las = LASSCF(mf, (4,),(4,), spin_sub=(1,))

# Localize the chosen fragment active spaces
#frag_atom_list = ((0,1),)
frag_atom_list = ((0,1),(2,3))
#frag_atom_list = ((0,1,2,3),)
loc_mo_coeff = las.localize_init_guess(frag_atom_list, mf.mo_coeff)

# Run LASSCF
las.kernel(loc_mo_coeff)
loc_mo_coeff = las.mo_coeff
print("LASSCF energy: ", las.e_tot)

ncore = las.ncore
ncas = las.ncas
ncas_sub = las.ncas_sub
nelec_cas = las.nelecas

# Checking subspace sizes and num. of electrons
print ("Ncore: ", ncore, "Ncas: ", ncas, "Ncas_sub: ", ncas_sub, "Nelec_cas: ", nelec_cas)

# CASCI h1 & h2 for VQE Hamiltonian
mc = mcscf.CASCI(mf,4,4)
mc.kernel(loc_mo_coeff)
cas_h1e, e_core = mc.h1e_for_cas()

eri_cas = mc.get_h2eff(loc_mo_coeff)
eri = ao2mo.restore(1, eri_cas,mc.ncas)

# Gets all acceptable operators for UCCSD
# excluding intra-fragment ones
def custom_excitations(num_spin_orbitals: int,
                           num_particles: Tuple[int, int],
                           num_sub: List[int]
                           ) -> List[Tuple[Tuple[int, ...], ...]]:
    excitations = []
    norb = int(num_spin_orbitals/2)
    uop = lasuccsd.gen_uccsd_op (norb, num_sub)
    a_idxs = uop.a_idxs
    i_idxs = uop.i_idxs
    for a,i in zip(a_idxs,i_idxs):
        excitations.append((tuple(i),tuple(a[::-1])))
    
    return excitations

hamiltonian = get_hamiltonian(None, mc.nelecas, mc.ncas, cas_h1e, eri)
#print(hamiltonian)

# Initialize using LASCI vector
# Code stolen from Riddhish
## This function makes a few assumptions
## 1. The civector is arranged as a 2D matrix of coeffs
##    of size [nalphastr, nbetastr]
## 2. The civector contains all configurations within
##    the (localized) active space
def get_so_ci_vec(ci_vec, nsporbs,nelec):
    lookup = {}
    cnt = 0
    norbs = nsporbs//2

    # Here, we set up a lookup dictionary which is
    # populated when either the number of alpha e-s
    # or the number of beta electrons is correct
    # It stores "bitstring" : decimal_value pairs
    ## The assumption is that nalpha==nbeta
    for ii in range (2**norbs):
        if f"{ii:0{norbs}b}".count('1') == np.sum(nelec)//2:
            lookup[f"{ii:0{norbs}b}"] = cnt
            cnt +=1
    # This is just indexing the hilber space from 0,1,...,mCn
    #print (lookup)

    # Here the spin orbital CI vector is populated
    # the same lookup is used for alpha and beta, but for two different
    # sections of the bitstring
    so_ci_vec = np.zeros(2**nsporbs)
    for kk in range (2**nsporbs):
        if f"{kk:0{nsporbs}b}"[norbs:].count('1')==nelec[0] and f"{kk:0{nsporbs}b}"[:norbs].count('1')==nelec[1]:
            so_ci_vec[kk] = ci_vec[lookup[f"{kk:0{nsporbs}b}"[norbs:]],lookup[f"{kk:0{nsporbs}b}"[:norbs]]]

    return so_ci_vec

qr1 = QuantumRegister(np.sum(ncas_sub)*2, 'q1')
new_circuit = QuantumCircuit(qr1)
new_circuit.initialize( get_so_ci_vec(las.ci[0][0],2*ncas_sub[0],las.nelecas_sub[0]) , [0,1,4,5])
new_circuit.initialize( get_so_ci_vec(las.ci[1][0],2*ncas_sub[1],las.nelecas_sub[1]) , [2,3,6,7])

# Gate counts for initialization
if args.dist == 0.0:
    target_basis = ['rx', 'ry', 'rz', 'h', 'cx']
    circ_for_counts = transpile(new_circuit, basis_gates=target_basis, optimization_level=0)
    init_op_dict = circ_for_counts.count_ops()
    init_ops = sum(init_op_dict.values())
    print("Operations: {}".format(init_op_dict))
    print("Total operations: {}".format(init_ops))

# Tracking the convergence of the VQE
counts = []
values = []
stds = []
def store_intermediate_result(eval_count, parameters, mean, std):
    counts.append(eval_count)
    values.append(mean)
    stds.append(std)

qubit_converter = QubitConverter(mapper = JordanWignerMapper(), two_qubit_reduction=False)
#init_test = HartreeFock(8,(2,2), qubit_converter)

# Setting up the VQE
backend = Aer.get_backend('aer_simulator')
#ansatz = custom_UCC(qubit_converter=qubit_converter, num_particles=(2,2), num_spin_orbitals=8, excitations=custom_excitations, initial_state=init_test, preserve_spin=False)
ansatz = custom_UCC(qubit_converter=qubit_converter, num_particles=(2,2), num_spin_orbitals=8, excitations=custom_excitations, initial_state=new_circuit, preserve_spin=False)
optimizer = L_BFGS_B(maxfun=10000, iprint=101)
init_pt = np.zeros(146)
#optimizer = COBYLA(maxiter=1000)
with Session(service = service , backend = backend):
    estimator = Estimator()
    options = Options()
    options.shots =  args.shots
    algorithm = VQE(ansatz=ansatz, optimizer=optimizer, quantum_instance=new_instance, initial_point=init_pt, callback=store_intermediate_result) 


# Gate counts for VQE (includes initialization)
if args.dist == 0.0:
    params = np.zeros(146)
    vqe_ops = 0
    circ_list = transpile(algorithm.construct_circuit(params, hamiltonian), basis_gates=target_basis, optimization_level=0)
    for circ in circ_list:
        vqe_op_dict = circ.count_ops()
        vqe_ops += sum(vqe_op_dict.values())
    print("Number of circuits in list: ",len(circ_list))
    print("Operations: {}".format(vqe_op_dict))
    print("Total operations: {}".format(vqe_ops))

# Running the VQE
t0 = time.time()
vqe_result = algorithm.compute_minimum_eigenvalue(hamiltonian)
print(vqe_result)
t1 = time.time()
print("Time taken for VQE: ",t1-t0)
print("VQE counts: ", counts)
print("VQE energies: ", values)

# Saving all relevant results in a dict
if args.dist == 0.0:
    np.save('results_{}_{}_vqe.npy'.format(args.shots, args.dist), {'init_op_dict': init_op_dict, 'init_ops': init_ops, 'vqe_op_dict':vqe_op_dict, 'vqe_ops': vqe_ops, 'vqe_result':vqe_result, 'vqe_en_vals':values, 'vqe_counts':counts, 'nuc_rep': las.energy_nuc()})
else:
    np.save('results_{}_{}_vqe.npy'.format(args.shots, args.dist), {'vqe_result':vqe_result, 'vqe_en_vals':values, 'vqe_counts':counts, 'nuc_rep': las.energy_nuc()})
