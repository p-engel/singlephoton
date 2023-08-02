# Routine for plotting the results of the g2 calculation


import numpy as np
import matplotlib.pyplot as plt


# handle data file
filenames = ['../data/g2_wpeak.txt', '../data/g2_wres.txt',
    '../data/g2_wdip.txt']
dataset = []
for filename in filenames:
    dataset.append(np.loadtxt(filename, skiprows=2))

tauc = 100854569.64740749   # [ps]
time = dataset[0][:,0]/tauc    # [units wrt tauc]

# plot
fig, axs = plt.subplots(4, sharex=True, figsize=(9, 9))
linestyles = [':','-.', '--']
labels = ['w=wpeak', 'w=res','w=dip']
axs[0].set(xlabel='t [tauc]', ylabel='g2(t, tau=0)')
# axs[0].set_ylim([0, 1.5])

for idx in range(len(labels)):
    axs[0].plot(time, dataset[idx][:,1], linestyles[idx], color='black', \
        lw=2, label=labels[idx])
    axs[0].legend()
    # Avg. photon number
    axs[idx+1].plot(time, dataset[idx][:,2], linestyles[idx], color='blue', \
        lw=2, label='Nc')
    axs[idx+1].plot(time, dataset[idx][:,3], linestyles[idx], color='red', 
        lw=2, label='Ne')
    axs[idx+1].set(xlabel='t [tauc]', ylabel='<n>(t, '+labels[idx]+' )')
    axs[idx+1].legend()

fig.suptitle('A discrete-cavity weakly coupled to' \
    + ' a weakly-driven broadband-emitter', fontsize=12)

plt.show()


# seventy-two and seventy-nine characters in a line for comment and code
# 34567890123456789012345678901234567890123456789012345678901234567890123456789
