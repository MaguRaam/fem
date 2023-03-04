# parse gmsh file and write the node.dat and ele.dat for fem solver:

import numpy as np
import gmshparser
import matplotlib.pylab as plt

mesh = gmshparser.parse("square.msh")

X = []
Y = []
T = []

# extract nodal coordinates:
for entity in mesh.get_node_entities():
    for node in entity.get_nodes():
        ncoords = node.get_coordinates()
        X.append(ncoords[0])
        Y.append(ncoords[1])

# extract connectivity of triangle elements:
for entity in mesh.get_element_entities():
    if entity.get_element_type() == 2: #triangle element 
        for element in entity.get_elements():
            conn = element.get_connectivity()
            T.append([conn[0] - 1, conn[1] - 1, conn[2] - 1])

# write node coordinates to file
np.savetxt("node.dat", np.column_stack([X, Y]), fmt="%.16f %.16f")

# write element connectivity to file
np.savetxt("ele.dat", T, fmt="%d %d %d")

# plot mesh
plt.figure()
plt.triplot(X, Y, T, '-r')
plt.axis('equal')
plt.tight_layout()
plt.savefig('square.png', dpi=300)
