import numpy as np
import matplotlib.pyplot as plt
from math import log

# log-log convergence plot of L2error and energy error:
def plot_convergence(h, error1, label1, error2, label2, filename):
    fig, ax = plt.subplots(figsize=(8,6))
    ax.loglog(h, error1, '-bo', label=label1, linewidth=2)
    ax.loglog(h, error2, '-ro', label=label2, linewidth=2)
    ax.set_xlabel(r'$h$', fontsize=16)
    ax.set_ylabel('Error', fontsize=16)
    ax.set_title('Convergence plot', fontsize=20)
    ax.set_xlim([min(h), max(h)])
    ax.legend(fontsize=14)
    fig.savefig(filename, dpi=300, bbox_inches='tight')

def convergence_rate(error, h):
    return [(log(error[i-1]) - log(error[i])) / (log(h[i - 1]) - log(h[i])) for i in range(1, np.size(error))]

# plot energy error of uh, and piu
def plot_energy_error(filename, save_as='energyerror_plot.png', dpi=300):
    h, energyerror_uh, energyerror_piu = np.loadtxt(filename, unpack=True)

    plt.plot(h, energyerror_uh, '-ro', label='$u_h$ ')
    plt.plot(h, energyerror_piu, '-go', label='$\pi u$ ')

    plt.title('Energy error', fontsize=16)
    plt.xlabel('$h$', fontsize=14)
    plt.ylabel('Energy error', fontsize=14)

    plt.legend(loc='best', fontsize=12)
    plt.style.use('seaborn-whitegrid')

    fig = plt.gcf()
    fig.set_size_inches(8, 6)
    fig.set_dpi(dpi)

    plt.rc('font', size=12)

    plt.savefig(save_as, dpi=dpi, bbox_inches='tight')


# plot energy functional:
import numpy as np
import matplotlib.pyplot as plt

def plot_energy_functional(data_path, plot_title, x_label, y_label):
    # load h and energy data
    h, energy_uexact, energy_uh, energy_piu = np.loadtxt(data_path, unpack=True)

    # plot h versus energy_uexact
    plt.plot(h, energy_uexact, '-bo', label='$u_{exact}$')

    # plot h versus energy_uh
    plt.plot(h, energy_uh, '-ro', label='$u_h$ ')

    # plot h versus energy_piu
    plt.plot(h, energy_piu, '-go', label='$\pi u$ ')

    # set plot title and labels
    plt.title(plot_title, fontsize=16)
    plt.xlabel(x_label, fontsize=14)
    plt.ylabel(y_label, fontsize=14)

    # set legend
    plt.legend(loc='best', fontsize=12)

    # set plot style
    plt.style.use('seaborn-whitegrid')

    # set figure size and dpi
    fig = plt.gcf()
    fig.set_size_inches(8, 6)
    fig.set_dpi(300)

    # set font size
    plt.rc('font', size=12)

    # save plot as high-resolution image
    plt.savefig('energy_plot.png', dpi=300, bbox_inches='tight')

    # show plot
    plt.show()

# load h and error 
h, l2error, energy_error = np.loadtxt("error.dat", unpack=True)

# plot convergence
plot_convergence(h, l2error, '$L_2$ error', energy_error, 'Energy error', "error.png")

# calculate and print convergence rates
l2_rates = convergence_rate(l2error, h)
energy_rates = convergence_rate(energy_error, h)
print("L2 convergence rates:", l2_rates)
print("Energy convergence rates:", energy_rates)
