import numpy as np
import matplotlib.pyplot as plt

x, Pf, f = np.loadtxt('sol.dat', unpack=True)

plt.plot(x, f, color='blue', marker = 'o', label='$sin(3\pi x)$')
plt.plot(x, Pf, color='red', marker = 'o', label='$P_h sin(3\pi x)$')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.legend()
plt.savefig('plot.png', dpi=500)
plt.clf()
