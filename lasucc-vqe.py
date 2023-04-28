Last login: Fri Apr 28 13:24:19 on ttys001
(base) omarelkashef@Omars-MacBook-Pro las-qpe % ssh omarelkashef@midway3.rcc.uchicago.edu
Password: 
Password: 
Duo two-factor login for omarelkashef

Enter a passcode or select one of the following options:

 1. Duo Push to +XX XXX XXX 0626
 2. Duo Push to Omar's IPad (iOS)
 3. SMS passcodes to +XX XXX XXX 0626

Passcode or option (1-3): 1
Success. Logging you in...
Activate the web console with: systemctl enable --now cockpit.socket

Last failed login: Fri Apr 28 16:23:36 CDT 2023 from 98.46.134.110 on ssh:notty
There was 1 failed login attempt since the last successful login.
Last login: Fri Apr 28 12:30:44 2023 from 10.150.123.135
[omarelkashef@midway3-login3 ~]$ squeue --me
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
[omarelkashef@midway3-login3 ~]$ cd /project
[omarelkashef@midway3-login3 project]$ cd lgagliardi/omarelkashef/
[omarelkashef@midway3-login3 omarelkashef]$ ls
9_500_03_10s   debugging-molcas.sh  new_v_2-5_7s_pdft  v_2-5_10t_1727940.err  v_2-5_10t_1727953.err  v_2-5_10t7s_cas       v_2-5_10t.err   v_2-5_10t_pt2      v_2-5_7s_pdft       vns3o_500_03_15s
9_500_03_15s   las-qpe              prep_calcs.sh      v_2-5_10t_1727940.out  v_2-5_10t_1727953.out  v_2-5_10t7s_cas.inp   v_2-5_10t.inp   v_2-5_10t.RICDLib  v_2-5_7s_pt2        vns3o_500_03_molc
9_500_03_molc  las-vqe              test               v_2-5_10t_1727946.err  v_2-5_10t_1727965.err  v_2-5_10t7s_HMC-PDFT  v_2-5_10t.out   v_2-5_10t.status   vanadium_align.xyz  xmldump
COPY           molcas.sh            v_2-5_10t          v_2-5_10t_1727946.out  v_2-5_10t_1727965.out  v_2-5_10t7s_pdft      v_2-5_10t_pdft  v_2-5_7s           vns3o_500_03_11s
[omarelkashef@midway3-login3 omarelkashef]$ cd las-qpe/
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe              aer_las_qpe_3834251.out  aer_las_qpe_3836130.out  fermionic_excitation_generator.py  las-qpe.py                qasm_las_qpe_3834749.out  transpiled_circuit.txt
aer_las_qpe_3834234.err  aer_las_qpe_3834678.err  aer_lasqpe.sh            get_geom.py                        __pycache__               qasm_las_qpe_3834751.err
aer_las_qpe_3834234.out  aer_las_qpe_3834678.out  custom_UCC.py            get_hamiltonian.py                 qasm_las_qpe              qasm_las_qpe_3834751.out
aer_las_qpe_3834251.err  aer_las_qpe_3836130.err  di-las-ucc.py            h4_sto3g_1.35296239.log            qasm_las_qpe_3834749.err  qasm_lasqpe.sh
[omarelkashef@midway3-login3 las-qpe]$ cd qasm_las_qpe/
[omarelkashef@midway3-login3 qasm_las_qpe]$ ls
[omarelkashef@midway3-login3 qasm_las_qpe]$ cd ..
[omarelkashef@midway3-login3 las-qpe]$ cd aer_las_qpe/
[omarelkashef@midway3-login3 aer_las_qpe]$ ls
[omarelkashef@midway3-login3 aer_las_qpe]$ cd ..
[omarelkashef@midway3-login3 las-qpe]$ vim aer_las_qpe_3836130.
[omarelkashef@midway3-login3 las-qpe]$ vim aer_las_qpe_3836130.err 
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe              aer_las_qpe_3834251.out  aer_las_qpe_3836130.out  fermionic_excitation_generator.py  las-qpe.py                qasm_las_qpe_3834749.out  transpiled_circuit.txt
aer_las_qpe_3834234.err  aer_las_qpe_3834678.err  aer_lasqpe.sh            get_geom.py                        __pycache__               qasm_las_qpe_3834751.err
aer_las_qpe_3834234.out  aer_las_qpe_3834678.out  custom_UCC.py            get_hamiltonian.py                 qasm_las_qpe              qasm_las_qpe_3834751.out
aer_las_qpe_3834251.err  aer_las_qpe_3836130.err  di-las-ucc.py            h4_sto3g_1.35296239.log            qasm_las_qpe_3834749.err  qasm_lasqpe.sh
[omarelkashef@midway3-login3 las-qpe]$ vim aer_lasqpe.sh 
[omarelkashef@midway3-login3 las-qpe]$ vim qasm_lasqpe.sh 
[omarelkashef@midway3-login3 las-qpe]$ vim qasm_lasqpe.sh 
[omarelkashef@midway3-login3 las-qpe]$ l
bash: l: command not found...
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe              aer_las_qpe_3834251.out  aer_las_qpe_3836130.out  fermionic_excitation_generator.py  las-qpe.py                qasm_las_qpe_3834749.out  transpiled_circuit.txt
aer_las_qpe_3834234.err  aer_las_qpe_3834678.err  aer_lasqpe.sh            get_geom.py                        __pycache__               qasm_las_qpe_3834751.err
aer_las_qpe_3834234.out  aer_las_qpe_3834678.out  custom_UCC.py            get_hamiltonian.py                 qasm_las_qpe              qasm_las_qpe_3834751.out
aer_las_qpe_3834251.err  aer_las_qpe_3836130.err  di-las-ucc.py            h4_sto3g_1.35296239.log            qasm_las_qpe_3834749.err  qasm_lasqpe.sh
[omarelkashef@midway3-login3 las-qpe]$ mkdir trash
[omarelkashef@midway3-login3 las-qpe]$ cp aer_las_qpe_3834234.err aer_las_qpe_3834234.out ./trash
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe              aer_las_qpe_3834251.out  aer_las_qpe_3836130.out  fermionic_excitation_generator.py  las-qpe.py                qasm_las_qpe_3834749.out  transpiled_circuit.txt
aer_las_qpe_3834234.err  aer_las_qpe_3834678.err  aer_lasqpe.sh            get_geom.py                        __pycache__               qasm_las_qpe_3834751.err  trash
aer_las_qpe_3834234.out  aer_las_qpe_3834678.out  custom_UCC.py            get_hamiltonian.py                 qasm_las_qpe              qasm_las_qpe_3834751.out
aer_las_qpe_3834251.err  aer_las_qpe_3836130.err  di-las-ucc.py            h4_sto3g_1.35296239.log            qasm_las_qpe_3834749.err  qasm_lasqpe.sh
[omarelkashef@midway3-login3 las-qpe]$ rm aer_las_qpe_3834234.err aer_las_qpe_3834234.out ./trash
rm: cannot remove './trash': Is a directory
[omarelkashef@midway3-login3 las-qpe]$ rm aer_las_qpe_3834234.err aer_las_qpe_3834234.out 
rm: cannot remove 'aer_las_qpe_3834234.err': No such file or directory
rm: cannot remove 'aer_las_qpe_3834234.out': No such file or directory
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe              aer_las_qpe_3834678.out  custom_UCC.py                      get_hamiltonian.py       qasm_las_qpe              qasm_las_qpe_3834751.out
aer_las_qpe_3834251.err  aer_las_qpe_3836130.err  di-las-ucc.py                      h4_sto3g_1.35296239.log  qasm_las_qpe_3834749.err  qasm_lasqpe.sh
aer_las_qpe_3834251.out  aer_las_qpe_3836130.out  fermionic_excitation_generator.py  las-qpe.py               qasm_las_qpe_3834749.out  transpiled_circuit.txt
aer_las_qpe_3834678.err  aer_lasqpe.sh            get_geom.py                        __pycache__              qasm_las_qpe_3834751.err  trash
[omarelkashef@midway3-login3 las-qpe]$ mv aer_las_qpe_3834251.err aer_las_qpe_3834251.out 
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe              aer_las_qpe_3834678.out  aer_lasqpe.sh  fermionic_excitation_generator.py  h4_sto3g_1.35296239.log  qasm_las_qpe              qasm_las_qpe_3834751.err  transpiled_circuit.txt
aer_las_qpe_3834251.out  aer_las_qpe_3836130.err  custom_UCC.py  get_geom.py                        las-qpe.py               qasm_las_qpe_3834749.err  qasm_las_qpe_3834751.out  trash
aer_las_qpe_3834678.err  aer_las_qpe_3836130.out  di-las-ucc.py  get_hamiltonian.py                 __pycache__              qasm_las_qpe_3834749.out  qasm_lasqpe.sh
[omarelkashef@midway3-login3 las-qpe]$ mv aer_las_qpe_3834251.out aer_las_qpe_3834251 ./trash
mv: cannot stat 'aer_las_qpe_3834251': No such file or directory
[omarelkashef@midway3-login3 las-qpe]$ mv aer_las_qpe_3834251.out ./trash
mv: cannot stat 'aer_las_qpe_3834251.out': No such file or directory
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe              aer_las_qpe_3836130.err  custom_UCC.py                      get_geom.py              las-qpe.py    qasm_las_qpe_3834749.err  qasm_las_qpe_3834751.out  trash
aer_las_qpe_3834678.err  aer_las_qpe_3836130.out  di-las-ucc.py                      get_hamiltonian.py       __pycache__   qasm_las_qpe_3834749.out  qasm_lasqpe.sh
aer_las_qpe_3834678.out  aer_lasqpe.sh            fermionic_excitation_generator.py  h4_sto3g_1.35296239.log  qasm_las_qpe  qasm_las_qpe_3834751.err  transpiled_circuit.txt
[omarelkashef@midway3-login3 las-qpe]$ mv aer_las_qpe_3834678.err aer_las_qpe_3834678.out ./trash
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe              aer_lasqpe.sh  fermionic_excitation_generator.py  h4_sto3g_1.35296239.log  qasm_las_qpe              qasm_las_qpe_3834751.err  transpiled_circuit.txt
aer_las_qpe_3836130.err  custom_UCC.py  get_geom.py                        las-qpe.py               qasm_las_qpe_3834749.err  qasm_las_qpe_3834751.out  trash
aer_las_qpe_3836130.out  di-las-ucc.py  get_hamiltonian.py                 __pycache__              qasm_las_qpe_3834749.out  qasm_lasqpe.sh
[omarelkashef@midway3-login3 las-qpe]$ mv aer_las_qpe_3836130.err aer_las_qpe_3836130.out ./trash
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe    di-las-ucc.py                      get_hamiltonian.py       __pycache__               qasm_las_qpe_3834749.out  qasm_lasqpe.sh
aer_lasqpe.sh  fermionic_excitation_generator.py  h4_sto3g_1.35296239.log  qasm_las_qpe              qasm_las_qpe_3834751.err  transpiled_circuit.txt
custom_UCC.py  get_geom.py                        las-qpe.py               qasm_las_qpe_3834749.err  qasm_las_qpe_3834751.out  trash
[omarelkashef@midway3-login3 las-qpe]$ mv qasm_las_qpe_3834749.err qasm_las_qpe_3834749.out ./trash
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe    custom_UCC.py  fermionic_excitation_generator.py  get_hamiltonian.py       las-qpe.py   qasm_las_qpe              qasm_las_qpe_3834751.out  transpiled_circuit.txt
aer_lasqpe.sh  di-las-ucc.py  get_geom.py                        h4_sto3g_1.35296239.log  __pycache__  qasm_las_qpe_3834751.err  qasm_lasqpe.sh            trash
[omarelkashef@midway3-login3 las-qpe]$ mv qasm_las_qpe_3834751.err qasm_las_qpe_3834751.out ./trash
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe    custom_UCC.py  fermionic_excitation_generator.py  get_hamiltonian.py       las-qpe.py   qasm_las_qpe    transpiled_circuit.txt
aer_lasqpe.sh  di-las-ucc.py  get_geom.py                        h4_sto3g_1.35296239.log  __pycache__  qasm_lasqpe.sh  trash
[omarelkashef@midway3-login3 las-qpe]$ cd trash
[omarelkashef@midway3-login3 trash]$ ls
aer_las_qpe_3834234.err  aer_las_qpe_3834251.out  aer_las_qpe_3834678.out  aer_las_qpe_3836130.out   qasm_las_qpe_3834749.out  qasm_las_qpe_3834751.out
aer_las_qpe_3834234.out  aer_las_qpe_3834678.err  aer_las_qpe_3836130.err  qasm_las_qpe_3834749.err  qasm_las_qpe_3834751.err
[omarelkashef@midway3-login3 trash]$ cd ..
[omarelkashef@midway3-login3 las-qpe]$ ls
aer_las_qpe    custom_UCC.py  fermionic_excitation_generator.py  get_hamiltonian.py       las-qpe.py   qasm_las_qpe    transpiled_circuit.txt
aer_lasqpe.sh  di-las-ucc.py  get_geom.py                        h4_sto3g_1.35296239.log  __pycache__  qasm_lasqpe.sh  trash
[omarelkashef@midway3-login3 las-qpe]$ vim di-las-ucc.py 

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
"di-las-ucc.py" 234L, 9425C written                                                                                                                                                    45,55          0%
