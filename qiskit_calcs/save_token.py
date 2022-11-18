from qiskit import IBMQ
from qiskit.providers.ibmq import least_busy

IBMQ.load_account() # Load account from disk
print(IBMQ.providers()  )  # List all available providers
provider = IBMQ.get_provider(hub='ibm-q')
print(provider.backends())

small_devices = provider.backends(filters=lambda x: x.configuration().n_qubits == 5
                                   and not x.configuration().simulator)
print(least_busy(small_devices))
