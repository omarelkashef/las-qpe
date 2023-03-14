# LAS-VQE

This code performs GB-LAS-UCC and DI-LAS-UCC calculations.

Requires:
1. Qiskit
2. Qiskit Nature (version < 0.5)
3. PySCF
4. MRH (github.com/MatthewRHermes/mrh)

GB-LAS-UCC steps:
1. Run `qpe.py` to obtain the ground state for each fragment, stored in `qpe_state.npy`
2. Run `las-vqe.py` using the initialize step, loading wave functions from `qpe_state.npy`

DI-LAS-UCC steps:
1. Run `las-vqe.py` using the initialize step and `get_so_ci_vector()`

HF-UCC steps:
1. Run `las-vqe.py` using the initialize step and the HartreeFock initial state
