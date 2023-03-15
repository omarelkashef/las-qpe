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

def get_hamiltonian(frag, nelecas_sub, ncas_sub, nuc_rep_en, h1_frag, h2_frag):
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
    print("Nuclear repulsion: ", nuc_rep_en)
    electronic_energy = ElectronicEnergy(
        [
            # Using MO basis here for simplified conversion
            OneBodyElectronicIntegrals(ElectronicBasis.MO, (h1_frag[frag], None)),
            TwoBodyElectronicIntegrals(ElectronicBasis.MO, (h2_frag[frag], h2_frag[frag], h2_frag[frag], None)),
        ],
        nuclear_repulsion_energy=nuc_rep_en,
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
    return hamiltonian

