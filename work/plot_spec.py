import matplotlib.pyplot as plt
import numpy as np

levels = np.loadtxt('levels.txt')
lengths = np.arange(4.0,10.01,.01)
for i, L in enumerate(lengths):
    plt.plot([L for a in levels[i,:]], levels[i,:],'.k',markersize=1)
plt.xlim([4,10])
plt.ylim([200,1500])
plt.ylabel(r'$E$ (MeV)')
plt.xlabel(r'L (fm)')
plt.show()
