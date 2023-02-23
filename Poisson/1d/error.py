import numpy as np
import matplotlib.pyplot as plt
from math import log

def plot_convergence(h, error, filename):
    fig, ax = plt.subplots()
    ax.loglog(h, error, '-bo', label='$L_{\infty}$ error')
    ax.set_xlabel(r'$h$')
    ax.set_ylabel(r'$E(h)$')
    ax.legend()
    fig.savefig(filename, dpi=300, bbox_inches='tight')


def convergence_rate(error):
    return [(log(error[i-1]) - log(error[i])) / log(2) for i in range(1, np.size(error))]

# load h and error abs(p - pexact)
h, error = np.loadtxt("error.dat", unpack=True)

# plot convergence
plot_convergence(h, error, "error.png")

# calculate and print convergence rates
rates = convergence_rate(error)
print("Convergence rates:", rates)
