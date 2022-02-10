import numpy as np
import matplotlib.pyplot as plt

an_list = [1,2,4,6,8,10]
frag_list = [0,1]

np_en = [-5.35375621, -5.35384458] 

#shots = 8192
shots = 10240

en_list = []
err_list = []
op_list = []

for frag in frag_list:
    for an in an_list:
        results = np.load('results_{}_{}.npy'.format(an, shots), allow_pickle=True).item()
        fig = plt.figure(figsize=(12,4.8))
        plt.bar(range(len(results['phases'][frag])), list(results['phases'][frag].values()))
        plt.xticks(range(len(results['phases'][frag])), list(results['phases'][frag].keys()), fontsize=8, rotation=45, ha="right")
        plt.xlabel('Ancilla qubit states')
        plt.ylabel('Probability')
        #fig.savefig("{}_qubits_{}_shots_c4h6_{}.png".format(an, shots, frag), bbox_inches="tight")

        most_likely_an = max(results['phases'][frag], key=results['phases'][frag].get)
        en_list.append(results['energies'][frag][most_likely_an])
        err_list.append(results['energies'][frag][most_likely_an] - np_en[frag])
        op_list.append(results['operations'][frag])

fig = plt.figure()
plt.plot(an_list, op_list[:6])
plt.xlabel('Number of ancilla qubits')
plt.ylabel('Number of gates')
#fig.savefig("gates_c4h6_0.png")

fig = plt.figure()
plt.plot(an_list, op_list[6:])
plt.xlabel('Number of ancilla qubits')
plt.ylabel('Number of gates')
#fig.savefig("gates_c4h6_1.png")

an_list = an_list[1:]
en_list = en_list[1:]
err_list = err_list[1:]

fig = plt.figure()
plt.plot(an_list, en_list[:5])
plt.axhline()
plt.xlabel('Number of ancilla qubits')
plt.ylabel('E [Hartree]')
fig.savefig("en_{}_shots_c4h6_0.png".format(shots), bbox_inches="tight")
fig = plt.figure()
plt.plot(an_list, err_list[:5])
plt.axhline()
plt.axhspan(0.0, 0.0016, facecolor='0.5', alpha=0.5)
plt.xlabel('Number of ancilla qubits')
plt.ylabel('E-E(FCI) [Hartree]')
fig.savefig("err_{}_shots_c4h6_0.png".format(shots), bbox_inches="tight")

