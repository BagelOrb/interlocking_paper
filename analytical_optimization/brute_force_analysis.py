import numpy as np

import matplotlib.pyplot as plt

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

wa_min = 2 * line_w
wb_min = 2 * line_w
w_max = 6 * wa_min
l_max = 6 * wa_min

h_min = 2 * layer_thickness
h_max = 6 * h_min

Nrz = 50
Nrw = 40
Nrx = 90
Nry = 90
shape = (Nrz, Nrw, Nrx, Nry)
# rzs = np.linspace(1.81, 1.82, Nrz)
# rws = np.linspace(1.09, 1.10, Nrw)
# rxs = np.linspace(1.80, 1.81, Nrx)
# rys = np.linspace(4.45, 4.50, Nry)
rzs = np.linspace(1, 5, Nrz)
rws = np.linspace(.5, 2, Nrw)
rxs = np.linspace(1, 5, Nrx)
rys = np.linspace(1, 5, Nry)

rz, rw, rx, ry = np.meshgrid(rzs, rws, rxs, rys, indexing='ij')

wa = np.full(shape, 2 * wa_min)
hc = np.full(shape, h_min)
wb = wa * rx
va = wa * rw
hf = hc * rz
vb = va * ry

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

viable_stress = stress
viable_stress[l > 3.6] = 0
best_idx = np.unravel_index(np.argmax(viable_stress), viable_stress.shape)

F_analytic = viable_stress * (wa + wb) * (hf + hc)
print(
    f"best: s: {viable_stress[best_idx]}\n"
    f"rz= {rz[best_idx]}, rw= {rw[best_idx]}, rx= {rx[best_idx]}, ry= {ry[best_idx]}, FF= {F_analytic[best_idx]}\n"
    f"wa={wa[best_idx]}, wb={wb[best_idx]}, va={va[best_idx]}, vb={vb[best_idx]}, hc={hc[best_idx]}, hf={hf[best_idx]}")

# fig, axs = plt.subplots(2, 2, figsize=(10, 6))
# axs[0, 0].scatter(l.reshape(-1), FEM_stress.reshape(-1), c=colors)
# axs[0, 0].set_title('Total length')
# plt.show()
