import numpy as np

import matplotlib.pyplot as plt
from matplotlib import pyplot
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

from math import sqrt

sa = 47
sb = 10.5
saz = sa # 33
sbz = sb # 9.0

ta = sa / sqrt(3)
tb = sb / sqrt(3)
taz = saz / sqrt(3)
tbz = sbz / sqrt(3)

line_w = .3
layer_thickness = .1

wa_min = line_w
wb_min = line_w
w_max = 6 * 2 * wa_min
l_max = 6 * 2 * wa_min

h_min = 2 * layer_thickness
h_max = 6 * h_min

# reading file

data_file = np.genfromtxt('IS-data-0527.csv', delimiter=',')
data_file[0, 0] = 0.2 # why doesn't it read that cell?!?!

Nhf = 5
Nlmax = 4
Nwb = Nva = 9
shape = (Nhf, Nlmax, Nwb, Nva)


FEM_stress = np.ndarray(shape=shape)
hf = np.ndarray(shape=shape)
lmax = np.ndarray(shape=shape)
wb = np.ndarray(shape=shape)
va = np.ndarray(shape=shape)


for hf_ in range(0, Nhf):
    hf[hf_, :, :, :] = data_file[(Nwb + 1) * hf_, 0]

for lmax_ in range(0, Nlmax):
    lmax[:, lmax_, :, :] = data_file[0, 1 + (Nva + 1) * lmax_]

for wb_ in range(0, Nwb):
    wb[:, :, wb_, :] = data_file[1 + wb_, 1]

for va_ in range(0, Nva):
    for lmax_ in range(0, Nlmax):
        va[:, lmax_:, :, va_] = data_file[0, 2 + lmax_ * (Nva + 1) +va_]


for hf_ in range(0, Nhf):
    for lmax_ in range(0, Nlmax):
        FEM_stress[hf_, lmax_, :, :] = data_file[1 + (Nwb + 1) * hf_:1 + (Nwb + 1) * hf_ + Nwb, 2 + (Nva + 1) * lmax_:2 + (Nva + 1) * lmax_ + Nva]

# FEM_stress[np.isnan(FEM_stress)] = 0
# print( np.any(np.isnan( FEM_stress)))

wa = np.full(shape, 2 * wa_min)
hc = np.full(shape, h_min)
vb = lmax - va
l = va + vb

F = FEM_stress * (wa + wb) * (hf + hc)

best_idx = np.unravel_index(np.argmax(FEM_stress), FEM_stress.shape)
print(f"best FEM: stess: {FEM_stress[best_idx]}")
print(f" wb={wb[best_idx]:.2f}; va={va[best_idx]:.2f}; lmax={lmax[best_idx]:.2f}; hf={hf[best_idx]:.2f}; ")
print(f" wa={wa[best_idx]:.2f}; vb={vb[best_idx]:.2f}; hc={hc[best_idx]:.2f}; F={F[best_idx]:.2f}")


bending = 0


gFs = []
gFs.append(sa * wa * hf)
gFs.append(sb * wb * hf)
gF_shear_a = 4 / 3 * hc * va * ta
gF_shear_b = 4 / 3 * hc * vb * tb
gF_bend_a = 2 * hc / wb * va * va * sa
gF_bend_b = 2 * hc / wa * vb * vb * sb
gF_cross_a = np.minimum(gF_shear_a, gF_bend_a)
gF_cross_b = np.minimum(gF_shear_b, gF_bend_b)
gFs.append(np.maximum(gF_cross_a, gF_cross_b))
# gFs.append(gF_cross_a)
# gFs.append(gF_cross_b)
gFs.append(taz * 4 / 3 * va * wa)
gFs.append(tbz * 4 / 3 * vb * wb)

minF = np.minimum.reduce(gFs)
stress = minF / ((wa + wb) * (hf + hc))
best_idx = np.unravel_index(np.argmax(stress), stress.shape)

F = stress * (wa + wb) * (hf + hc)
print(f"best ana: stress: {stress[best_idx]}")
print(f" wb={wb[best_idx]:.2f}; va={va[best_idx]:.2f}; lmax={l_max:.2f}; hf={hf[best_idx]:.2f}; ")
print(f" wa={wa[best_idx]:.2f}; vb={vb[best_idx]:.2f}; hc={hc[best_idx]:.2f}; F={F[best_idx]:.2f}")

names = ['wa', 'wb', 'va', 'vb', 'hf', 'hc', 'design', 'tension a', 'tension b', 'cross', 'shear Z a', 'shear Z b']
cross_gs_names = ['shear a', 'shear b', 'bend a', 'bend b']

print("\n== prediction ratio per failure mode ==")
prediction_ratio = stress / FEM_stress
for i, gF in enumerate(gFs):
    ratios = prediction_ratio[minF == gF]
    if ratios.size > 0:
        print(f"g{i} {names[i+7]}:  {np.average(ratios):.3f}, stdev: {np.std(ratios):.3f}")
    else:
        print(f"g{i} {names[i+7]}")

print("- cross beam sub failure modes -")
print(f"{np.average(prediction_ratio[minF == gF_shear_a]):.3f}, stdev: {np.std(prediction_ratio[minF == gF_shear_a]):.3f}")
print(f"{np.average(prediction_ratio[minF == gF_shear_b]):.3f}, stdev: {np.std(prediction_ratio[minF == gF_shear_b]):.3f}")
print(f"{np.average(prediction_ratio[minF == gF_bend_a]) :.3f}, stdev: {np.std(prediction_ratio[minF == gF_bend_a]):.3f}")
print(f"{np.average(prediction_ratio[minF == gF_bend_b]) :.3f}, stdev: {np.std(prediction_ratio[minF == gF_bend_b]):.3f}")




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
    F / hc * np.minimum(np.maximum(3/4 * sqrt(3), wb / ( 2 * va)) / (va * sa), np.maximum(3/4 * sqrt(3), wa / ( 2 * vb)) / (vb * sb)) - 1)  # cross
gs.append(1 * F / (4 * va * wa * taz) - 1)  # z shear
gs.append(1 * F / (4 * vb * wb * tbz) - 1)

cross_gs = [
    3 * F / (4 * va * hc * ta) - 1  # shear
    , 3 * F / (4 * vb * hc * tb) - 1
    , F * wb / (2 * va * va * hc * sa) - 1  # bending
    , F * wa / (2 * vb * vb * hc * sb) - 1]


print("\n== Active constraints: ", end="")
for i, g in enumerate(gs):
    if abs(g[best_idx]) < .00001:
        print(names[i], end=",")
print("")

print("-- cross beam constraints: ", end="")
for i, g in enumerate(cross_gs):
    if abs(g[best_idx]) < .00001:
        print(cross_gs_names[i], end=",")
print("\n")

constraints = np.maximum.reduce(gs)
for i, g in enumerate(gs):
    if g.max() > .001:
        print(f"constraint g{i} is violated! : {g.max():.4f}")



colors = np.zeros((Nhf * Nlmax * Nva * Nwb, 3))
colors[:, 0] = lmax.reshape(-1) / lmax.max()
colors[:, 1] = wb.reshape(-1) / wb.max()
colors[:, 2] = hf.reshape(-1) / hf.max()
# colors[:, 2] = va.reshape(-1) / va.max()

fig, axs = plt.subplots(2, 2, figsize=(10, 6))

axs[0, 0].scatter(l.reshape(-1), FEM_stress.reshape(-1), c=colors)
axs[0, 0].set_title('Total length')
axs[0, 1].scatter(hf.reshape(-1), FEM_stress.reshape(-1), c=colors)
axs[0, 1].set_title('hf')
axs[1, 0].scatter(wb.reshape(-1), FEM_stress.reshape(-1), c=colors)
axs[1, 0].set_title('wb')
# axs[1, 1].scatter(va.reshape(-1), FEM_stress.reshape(-1), c=colors)
# axs[1, 1].set_title('va')
axs[1, 1].scatter(lmax.reshape(-1), FEM_stress.reshape(-1), c=colors)
axs[1, 1].set_title('lmax')

for ax in axs.flat:
    ax.set(ylabel='max stress')
for ax in axs.flat:
    ax.label_outer()

# plt.show()

#

idx = best_idx # (2, 2, 2, 2)

fig, axs = plt.subplots(2, 3, figsize=(10, 6), subplot_kw={'projection': '3d'})
# axs[0, 0].plot_trisurf(hf[:, :, idx[2], idx[3]].reshape(-1), wb[:, :, idx[2], idx[3]].reshape(-1), FEM_stress[:, :, idx[2], idx[3]].reshape(-1), edgecolor= 'none')
axs[0, 1].plot_surface(hf[:, idx[1], :, idx[3]], wb[:, idx[1], :, idx[3]], FEM_stress[:, idx[1], :, idx[3]], edgecolor='none')
axs[0, 2].plot_surface(hf[:, idx[1], idx[2], :], wb[:, idx[1], idx[2], :], FEM_stress[:, idx[1], idx[2], :], edgecolor='none')
axs[1, 0].plot_surface(hf[idx[0], :, :, idx[3]], wb[idx[0], :, :, idx[3]], FEM_stress[idx[0], :, :, idx[3]], edgecolor='none')
axs[1, 1].plot_surface(hf[idx[0], :, idx[2], :], wb[idx[0], :, idx[2], :], FEM_stress[idx[0], :, idx[2], :], edgecolor='none')
# axs[1, 2].plot_surface(hf[idx[0], idx[1], :, :], wb[idx[0], idx[1], :, :], FEM_stress[idx[0], idx[1], :, :], edgecolor='none')

# plt.show()

# print(np.maximum.reduce(gs))
