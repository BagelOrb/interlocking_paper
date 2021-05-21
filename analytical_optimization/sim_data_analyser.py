import numpy as np

import matplotlib.pyplot as plt
from matplotlib import pyplot

data_file = np.genfromtxt('../../../simulation/IS90-data.csv', delimiter=',')

sa = 47
sb = 10.5
saz = 33
sbz = 9.0

ta = .5 * sa
tb = .5 * sb
taz = .5 * saz
tbz = .5 * sbz

line_w = .3
layer_thickness = .1

wa_min = line_w
wb_min = line_w
w_max = 6 * 2 * wa_min
l_max = 6 * 2 * wa_min

h_min = 2 * layer_thickness
h_max = 6 * h_min

Nrz = 5
Nrw = 4
Nrx = Nry = 9
shape = (Nrz, Nrw, Nrx, Nry)

FEM_stress = np.ndarray(shape=shape)
rz = np.ndarray(shape=shape)
rx = np.ndarray(shape=shape)
ry = np.ndarray(shape=shape)
rw = np.ndarray(shape=shape)

for rz_ in range(0, Nrz):
    rz[rz_, :, :, :] = data_file[0, 10 * rz_]

for rw_ in range(0, Nrw):
    rw[:, rw_, :, :] = data_file[1 + 10 * rw_, 0]

for rx_ in range(0, Nrx):
    rx[:, :, rx_, :] = data_file[2 + rx_, 0]

for ry_ in range(0, Nry):
    ry[:, :, :, ry_] = data_file[1, 1 + ry_]

for rz_ in range(0, 5):
    for rw_ in range(0, 4):
        FEM_stress[rz_, rw_, 0:9, 0:9] = data_file[2 + 10 * rw_:2 + 10 * rw_ + 9, 1 + 10 * rz_:1 + 10 * rz_ + 9]
FEM_stress[np.isnan(FEM_stress)] = 0

def toCsv(data, filename):
    output_data = data_file
    for rz_ in range(0, 5):
        for rw_ in range(0, 4):
            output_data[2 + 10 * rw_:2 + 10 * rw_ + 9, 1 + 10 * rz_:1 + 10 * rz_ + 9] = data[rz_, rw_, 0:9, 0:9]
    np.savetxt(filename, output_data, delimiter=',')


wa = np.full(shape, 2 * wa_min)
hc = np.full(shape, h_min)
wb = wa * rx
va = wa * rw
hf = hc * rz
vb = va * ry


best_idx = np.unravel_index(np.argmax(FEM_stress), FEM_stress.shape)
print(f"best: stess={FEM_stress[best_idx]}")
print(f" rx={rx[best_idx]:.4f}; ry={ry[best_idx]:.4f}; rw={rw[best_idx]:.4f}; rz={rz[best_idx]:.4f}; ")
print(f" wa={wa[best_idx]:.4f}; wb={wb[best_idx]:.4f}; va={va[best_idx]:.4f}; vb={vb[best_idx]:.4f}; hf={hf[best_idx]:.4f}; hc={hc[best_idx]:.4f}; ")

F = FEM_stress * (wa + wb) * (hf + hc)

gs = []
gs.append(1 - wa / 2 / wa_min)
gs.append(1 - wb / 2 / wb_min)
gs.append(1 - va / wa_min)
gs.append(1 - vb / wb_min)
gs.append(1 - hf / h_min)
gs.append(1 - hc / h_min)
gs.append((va + vb) / l_max - 1)
gs.append(F / (wa * hf * sa) - 1)  # tensile
gs.append(F / (wb * hf * sb) - 1)
gs.append(
    3 * F / (4 * hc) * np.minimum(np.maximum(2, wb / va) / (va * sa), np.maximum(2, wa / vb) / (vb * sb)) - 1)  # cross
gs.append(3 * F / (4 * va * wa * taz) - 1)  # z shear
gs.append(3 * F / (4 * vb * wa * tbz) - 1)

l = va + vb


constraints = np.maximum.reduce(gs)
toCsv(constraints, 'constraints.csv')

mech_constraints = np.maximum.reduce(gs[6:14])
toCsv(mech_constraints, 'mech_constraints.csv')


colors = np.zeros((Nrz * Nrw * Nry * Nrx, 3))
colors[:, 0] = rw.reshape(-1) / rw.max()
colors[:, 1] = rx.reshape(-1) / rx.max()
colors[:, 2] = rz.reshape(-1) / rz.max()
# colors[:, 2] = ry.reshape(-1) / ry.max()

fig, axs = plt.subplots(2, 2, figsize=(10, 6))

axs[0, 0].scatter(l.reshape(-1), FEM_stress.reshape(-1), c=colors)
axs[0, 0].set_title('Total length')
axs[0, 1].scatter(rz.reshape(-1), FEM_stress.reshape(-1), c=colors)
axs[0, 1].set_title('rz')
axs[1, 0].scatter(rx.reshape(-1), FEM_stress.reshape(-1), c=colors)
axs[1, 0].set_title('rx')
# axs[1, 1].scatter(ry.reshape(-1), FEM_stress.reshape(-1), c=colors)
# axs[1, 1].set_title('ry')
axs[1, 1].scatter(rw.reshape(-1), FEM_stress.reshape(-1), c=colors)
axs[1, 1].set_title('rw')

for ax in axs.flat:
    ax.set(ylabel='max stress')
for ax in axs.flat:
    ax.label_outer()

#plt.show()

# print(np.maximum.reduce(gs))
