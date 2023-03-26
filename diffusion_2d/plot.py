import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

# Read data from file
with open('sol.dat', 'r') as f:
    nnodes, nelems = map(int, f.readline().split())

    coords = np.zeros((nnodes, 2))
    for i in range(nnodes):
        coords[i] = np.array(list(map(float, f.readline().split())))

    conn = np.zeros((nelems, 3), dtype=int)
    for i in range(nelems):
        conn[i] = np.array(list(map(int, f.readline().split())))

    sol = np.array([float(f.readline()) for i in range(nnodes)])

# Construct triangulation
tri = Triangulation(coords[:, 0], coords[:, 1], conn)

# Plot solution
fig, ax = plt.subplots(figsize=(6, 4))
tcf = ax.tripcolor(tri, sol, cmap='plasma', edgecolors='none', linewidths=0.2)

# Add colorbar
cbar = fig.colorbar(tcf, ax=ax)
cbar.ax.set_ylabel('', rotation=270, labelpad=10)

# Set axis labels and limits
ax.set(xlabel='X', ylabel='Y', xlim=(coords[:, 0].min(), coords[:, 0].max()), 
       ylim=(coords[:, 1].min(), coords[:, 1].max()), title='Solution')

# Save plot to file
plt.savefig('solution.png', dpi=300, bbox_inches='tight')
