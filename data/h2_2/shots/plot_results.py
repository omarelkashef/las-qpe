import numpy as np
import matplotlib.pyplot as plt

an_list = [1,2,4,6,8,10]
frag_list = [0,1]
shots_list = [1024, 2048, 4096, 6144, 8192, 10240]
an = 10

np_en = [-1.62957716, -1.62941649] 

#shots = 8192

en_list = []
err_list = []
max_prob_list = []
op_list = []

for frag in frag_list:
    for shots in shots_list:
        results = np.load('results_{}_{}.npy'.format(an, shots), allow_pickle=True).item()
        #fig = plt.figure(figsize=(12,4.8))
        #plt.bar(range(len(results['phases'][frag])), list(results['phases'][frag].values()))
        #plt.xticks(range(len(results['phases'][frag])), list(results['phases'][frag].keys()), fontsize=8, rotation=45, ha="right")
        #plt.xlabel('Ancilla qubit states')
        #plt.ylabel('Probability')
        #fig.savefig("{}_qubits_{}_shots_h2_{}.png".format(an, shots, frag), bbox_inches="tight")

        most_likely_an = max(results['phases'][frag], key=results['phases'][frag].get)
        max_prob_list.append(max(results['phases'][frag].values()))
        en_list.append(results['energies'][frag][most_likely_an])
        err_list.append(results['energies'][frag][most_likely_an] - np_en[frag])
        op_list.append(results['operations'][frag])

fig = plt.figure()
plt.plot(shots_list, max_prob_list[:6])
plt.xlabel('Number of measurements')
plt.ylabel('Probability of max probability state')
fig.savefig("prob_shots_h2_0.png")

fig = plt.figure()
plt.plot(shots_list, max_prob_list[6:])
plt.xlabel('Number of measurements')
plt.ylabel('Probability of max probability state')
fig.savefig("prob_shots_h2_1.png")

fig = plt.figure()
plt.plot(shots_list, op_list[:6])
plt.xlabel('Number of measurements')
plt.ylabel('Number of gates')
fig.savefig("gates_shots_h2_0.png")

fig = plt.figure()
plt.plot(shots_list, op_list[6:])
plt.xlabel('Number of measurements')
plt.ylabel('Number of gates')
fig.savefig("gates_shots_h2_1.png")
