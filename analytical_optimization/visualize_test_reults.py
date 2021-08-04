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

plot = ax.scatter(wb, va, hf, marker='o', edgecolors='black', c=stress, cmap=plt.get_cmap('jet'))

fig.colorbar(plot, fraction=0.046*.5, pad=0.1)

ax.set_xlabel('$w_b$')
ax.set_ylabel('$v_a$')
ax.set_zlabel('$h_f$')

#ax.set_xticks(np.arange(wb.min(), wb.max(), .3))
ax.set_xticks(np.unique(wb))
ax.set_yticks(np.arange(va.min(), va.max(), .3), minor=True)
ax.set_yticks(np.arange(va.min(), va.max(), .6), minor=False)
#ax.set_yticks(np.unique(va))
#ax.set_zticks(np.unique(hf))

plt.margins(0)

plt.show()