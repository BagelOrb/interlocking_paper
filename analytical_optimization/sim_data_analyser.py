import numpy as np

import matplotlib.pyplot as plt
from matplotlib import pyplot

data_file = np.genfromtxt('../../simulation/IS90-data.csv', delimiter=',')

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
gs.append(3 * F / (4 * va * hc * ta) - 1)  # shear
gs.append(3 * F / (4 * vb * hc * tb) - 1)
gs.append(3 * F / (4 * va * wa * taz) - 1)  # z shear
gs.append(3 * F / (4 * vb * wa * tbz) - 1)
gs.append(3 * F * wb / (4 * va * va * hc * sa) - 1)  # bending
gs.append(3 * F * wa / (4 * vb * vb * hc * sb) - 1)

gFs = []
gFs.append(sa * wa * hf)
gFs.append(sb * wb * hf)
gFs.append(.5 * sa * 4 / 3 * va * hc)
gFs.append(.5 * sb * 4 / 3 * vb * hc)
gFs.append(taz * 4 / 3 * va * wa)
gFs.append(tbz * 4 / 3 * vb * wb)
gFs.append(sa * 4 / 3 * va * va * hc / wb)
gFs.append(sb * 4 / 3 * vb * vb * hc / wa)

l = va + vb

minF = np.minimum.reduce(gFs)
stress = minF / ((wa + wb) * (hf + hc))
toCsv(stress, 'stress.csv')

viable_stress = stress
viable_stress[l > 3.6] = 0
F_analytic = viable_stress * (wa + wb) * (hf + hc)
best_idx = np.unravel_index(np.argmax(viable_stress), viable_stress.shape)
print(f"best: s: {viable_stress[best_idx]}, rz: {rz[best_idx]}, rw: {rw[best_idx]}, rx: {rx[best_idx]}, ry: {ry[best_idx]}, F: {F_analytic[best_idx]}")


constraints = np.maximum.reduce(gs)
toCsv(constraints, 'constraints.csv')

mech_constraints = np.maximum.reduce(gs[7:15])
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
