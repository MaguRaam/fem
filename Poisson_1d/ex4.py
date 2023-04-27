#1d fem code to solve vertical deflection of a bar:

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve


# material parameters
Ac = 0.1                # cross sectional area (m2)
E = 200.0e9             # young's modulus  (N/m2)
P = 100.0               # force per unit length (N/m)

# local stiffness matrix:
def klocal(a, b):
    h = b - a
    return np.array([[1.0/h, -1.0/h], [-1.0/h, 1.0/h]])

# compute local force vector:
def flocal(a, b):
    h = b - a
    return np.array([0.5 * h, 0.5 * h])

def exact_solution(x):
        return (0.5*P*(x - 10.0)*x)/(Ac*E)

# create grid points to define Vh space on [0, 10m] domain:
L = 10.0
nElements = 20
nNodes = nElements + 1
dh = L/nElements
x = np.linspace(0, L, nNodes)

# initialize global stiffness matrix and force vector
K = lil_matrix((nNodes, nNodes))
F = np.zeros(nNodes)

# assemble the global stiffness matrix and force vector
for e in range(nElements):
    
    a, b = x[e], x[e + 1]       # element endpoints
    ke = -Ac*E*klocal(a, b)     # local stiffness matrix
    fe = P*flocal(a, b)         # local force vector

    # add the element contributions to the global matrices
    for i in range(2):
        for j in range(2):
            K[e + i, e + j] += ke[i, j]

    F[e:e+2] += fe

# apply homogeneous dirichlet boundary conditions

# fixed left endpoint
K[0, :] = 0
K[0, 0] = 1
F[0] = 0

# fixed right endpoint
K[-1, :] = 0
K[-1, -1] = 1
F[-1] = 0

# convert the stiffness matrix to csr format for efficient solving
K = K.tocsr()

# solve the system of equations for the nodal displacements
u = spsolve(K, F)
uexact = exact_solution(x)

# post processing:

# create the plot
plt.plot(x, u, '-r', label='FEM solution')
plt.plot(x, uexact, 'o', label='Exact solution')

# set axis labels and title
plt.xlabel('Position (m)')
plt.ylabel('Deflection (m)')
plt.title('Vertical Deflection of a Bar')
plt.legend()

# save the plot in a high-quality image format (e.g., PNG or PDF)
plt.savefig('deflection_plot.png', dpi=300, bbox_inches='tight')
