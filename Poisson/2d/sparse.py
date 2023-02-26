import scipy.io as sio
import matplotlib.pyplot as plt


def sparsity_pattern(filename):
    # Load the matrix from the Matrix Market file
    mat = sio.mmread(filename)

    # Plot the matrix
    plt.spy(mat)

    # Save the plot to a high-quality image file
    plt.savefig('K.png', dpi=300)


sparsity_pattern("K.mtx")