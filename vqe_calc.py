import numpy as np
from qiskit.algorithms.minimum_eigensolvers import VQE , NumPyMinimumEigensolver
from qiskit.algorithms.optimizers import SLSQP
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.circuit.library import UCCSD , HartreeFock
from qiskit_nature.second_q.mappers import ParityMapper , QubitConverter
from qiskit.primitives import Estimator
from qiskit_nature.second_q.transformers import FreezeCoreTransformer


# Define molecule
molecule = 'H .0 .0 -{0}; H .0 .0 {0}'

distances = np.arange(0.01,2.5,0.3)
vqe_energies = []
hf_energies = []
exact_energies = []

for i , distance in enumerate(distances):
    print("step {}".format(i))
    
    #set up experiment 
    driver = PySCFDriver(atom = molecule.format(distance/2))
    transformer = FreezeCoreTransformer()
    q_molecule = transformer.transform(driver.run())
    num_particles = q_molecule.num_particles
    num_spatial_orbitals = q_molecule.num_spatial_orbitals
    qubit_converter = QubitConverter(ParityMapper(), two_qubit_reduction=True)
    qubit_hamiltonian = qubit_converter.convert(q_molecule.second_q_ops()[0],
                                               num_particles = num_particles)



    #exact eigensolver
    exact_solver = NumPyMinimumEigensolver().compute_minimum_eigenvalue(qubit_hamiltonian)
    exact_energy = q_molecule.interpret(exact_solver).groundenergy
    exact_energies.append(exact_energy)

    #VQE
    estimator = Estimator()
    initial_state = HartreeFock(num_spatial_orbitals, num_particles,
                                qubit_converter)
    hf_energy = estimator.run([initial_state], [qubit_hamiltonian],
                              []).result().values[0]
    hf_energies.append(hf_energy)
    ansatz = UCCSD(num_spatial_orbitals , num_particles , qubit_converter,
                   initial_state=initial_state) 
    optimizer = SLSQP(maxiter=1000)
    vqe = VQE(estimator, ansatz, optimizer)
    vqe.initial_point = np.zeros(ansatz.num_parameters)
    vqe_solver = vqe.compute_minimum_eigenvalue(qubit_hamiltonian)
    vqe_energy = q_molecule.interpret(vqe_solver).groundenergy
    vqe_energies.append(vqe_energy)
    print(hf_energy,vqe_energy,exact_energy)