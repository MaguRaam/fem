import matplotlib.pyplot as plt
import numpy as np
import triangle as tr

# vertices:
v = [[0, 0], [0, 1], [1, 1], [1, 0]]

# triangulate:
t = tr.triangulate({'vertices': v}, 'qa0.001')

boundary_nodes = set()
for triangle in t['triangles']:
    for vertex in triangle:
        if t['vertex_markers'][vertex] == 1:
            boundary_nodes.add(vertex)

# Write the coordinates to a file
with open('coordinates.dat', 'w') as f:
    f.write('{}\n'.format(len(t['vertices'])))
    for i, vertex in enumerate(t['vertices']):
        f.write('{} {}\n'.format(vertex[0], vertex[1]))

# Write the connectivity to a file
with open('connectivity.dat', 'w') as f:
    f.write('{}\n'.format(len(t['triangles'])))
    for triangle in t['triangles']:
        f.write('{} {} {}\n'.format(triangle[0], triangle[1], triangle[2]))

# Write the unique boundary node numbers to a file
with open('boundarynodes.dat', 'w') as f:
    f.write('{}\n'.format(len(boundary_nodes)))
    for node in boundary_nodes:
        f.write('{}\n'.format(node))

# plot the mesh:
plt.figure()
plt.triplot(t['vertices'][:, 0], t['vertices'][:, 1], t['triangles'], '-k')
plt.scatter(t['vertices'][:, 0], t['vertices'][:, 1], s=1, c = 'black')
plt.axis('equal')
plt.tight_layout()
plt.savefig('mesh.png', dpi=300)