import numpy as np
import matplotlib.pyplot as plt

an_list = [1,2,4,6,8,10]
frag_list = [0,1]

np_en = [-1.62957716, -1.62941649] 

#shots = 8192
shots = 10240

en_list = []
err_list = []
op_list = []

for frag in frag_list:
    for an in an_list:
        results = np.load('results_{}_{}.npy'.format(an, shots), allow_pickle=True).item()
        fig = plt.figure(figsize=(12,4.8))
        state_list = list(results['phases'][frag].keys())
        prob_list = list(results['phases'][frag].values())
        new_state_list = []
        new_prob_list = []
        for idx, prob in enumerate(prob_list):
            if prob < 1e-4:
                print("Prob {} too small".format(prob))
            else:
                new_prob_list.append(prob)
                new_state_list.append(state_list[idx])
        print(prob_list,new_prob_list)
        print(state_list,new_state_list)
        plt.bar(range(len(new_state_list)), new_prob_list)
        plt.xticks(range(len(new_state_list)), new_state_list, fontsize=12, rotation=45, ha="right")
        plt.xlabel('Ancilla qubit states', fontsize=16)
        plt.ylabel('Probability', fontsize=16)
        fig.savefig("{}_qubits_{}_shots_h2_{}.png".format(an, shots, frag), bbox_inches="tight", dpi=300)

        most_likely_an = max(results['phases'][frag], key=results['phases'][frag].get)
        en_list.append(results['energies'][frag][most_likely_an])
        err_list.append(results['energies'][frag][most_likely_an] - np_en[frag])
        op_list.append(results['operations'][frag])

fig = plt.figure()
plt.plot(an_list, op_list[:6])
plt.xlabel('Number of ancilla qubits')
plt.ylabel('Number of gates')
fig.savefig("gates_h2_0.png")

fig = plt.figure()
plt.plot(an_list, op_list[6:])
plt.xlabel('Number of ancilla qubits')
plt.ylabel('Number of gates')
fig.savefig("gates_h2_1.png")

#an_list = an_list[1:]
#en_list = en_list[1:]
#err_list = err_list[1:]

fig = plt.figure()
plt.plot(an_list, en_list[6:], 'o-')
plt.axhline(color='k')
plt.xlabel('Number of ancilla qubits', fontsize=16)
plt.ylabel('E [Hartree]', fontsize=16)
plt.ylim(-0.4, 0.8)
fig.savefig("en_{}_shots_h2_1.png".format(shots), bbox_inches="tight", dpi=300)
fig = plt.figure()
plt.plot(an_list, err_list[6:], 'o-')
plt.axhline(color='k')
plt.axhspan(-0.355, 0.355, facecolor='0.5', alpha=0.5)
plt.xlabel('Number of ancilla qubits', fontsize=16)
plt.ylabel('E-E(FCI) [Hartree]', fontsize=16)
plt.ylim(-0.4, 0.8)
fig.savefig("err_{}_shots_h2_1.png".format(shots), bbox_inches="tight", dpi=300)

