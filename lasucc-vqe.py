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

parser = ArgumentParser(description='Do LAS-VQE, specifying num of ancillas and shots')
parser.add_argument('--dist', type=float, default=1.35296239, help='distance of H2s from one another or C=C bond distance scaling in butadiene')
parser.add_argument('--shots', type=int, default=1024, help='number of shots for the simulator')
args = parser.parse_args()

# Prints info at every VQE iteration
logging.basicConfig(level='DEBUG')

# Define molecule: (H2)_2
xyz = get_geom('far', dist=args.dist)
#xyz = '''H 0.0 0.0 0.0
#             H 1.0 0.0 0.0
#             H 0.2 1.6 0.1
"di-las-ucc.py" 234L, 9425C written                                                                                                                                                    23,33          0%
