import numpy as np
import matplotlib.pyplot as plt

x, u, uexact = np.loadtxt('sol.dat', unpack=True)

plt.plot(x, u, 'bo',label='numerical sol.')
plt.plot(x, uexact, '-r', label='exact sol.')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.legend()
plt.savefig('plot.png', dpi=500)
plt.clf()