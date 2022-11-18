from qiskit import IBMQ
from qiskit.providers.ibmq import least_busy

#IBMQ.save_account('f101fc398a5b97b9dda60c4528bb284a46f8cd2a3a8cfbfcd5f55866f4bd0161f94c15a015f27e38c17c9aa65d7369be5d48e2544717420c21718653b7411cb6')
IBMQ.load_account() # Load account from disk
print(IBMQ.providers()  )  # List all available providers
provider = IBMQ.get_provider(hub='ibm-q')
print(provider.backends())

small_devices = provider.backends(filters=lambda x: x.configuration().n_qubits == 5
                                   and not x.configuration().simulator)
print(least_busy(small_devices))
