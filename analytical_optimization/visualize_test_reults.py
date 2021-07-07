import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



data_file = np.genfromtxt('../../../testing/results.csv', delimiter=';')

wb = data_file[2:,2]
va = data_file[2:,3]
lmax = data_file[2:,4]
hf = data_file[2:,5]
stress = data_file[2:,13]  # divided by 5x5 cells
stress2 = data_file[2:,15]  # divided by total cross section area


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

plot = ax.scatter(wb, va, hf, marker='o', edgecolors='black', c=stress, cmap=plt.get_cmap('rainbow'))

fig.colorbar(plot, fraction=0.046*.5, pad=0.1)

ax.set_xlabel('$w_b$')
ax.set_ylabel('$v_a$')
ax.set_zlabel('$h_f$')

plt.show()