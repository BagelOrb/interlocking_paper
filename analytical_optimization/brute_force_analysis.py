import numpy as np

import matplotlib.pyplot as plt
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

N = 100
Nwb = N
Nva = N
Nhf = N
shape = (Nwb, Nva, Nhf)
wbs = np.linspace(2 * wb_min, w_max - 2 * wa_min, Nwb)
vas = np.linspace(wa_min, l_max - wb_min, Nva)
hfs = np.linspace(h_min, h_max - h_min, Nhf)

if True:
    wb_opt = 1.98
    va_opt = 2.07
    hf_opt = 0.72
    d = .1
    wbs = np.linspace(wb_opt - d, wb_opt + d, Nwb)
    vas = np.linspace(va_opt - d, va_opt + d, Nva)
    hfs = np.linspace(hf_opt - d, hf_opt + d, Nhf)

wb, va, hf = np.meshgrid(wbs, vas, hfs, indexing='ij')

l = np.full(shape, l_max)
wa = np.full(shape, 2 * wa_min)
hc = np.full(shape, h_min)
vb = l - va

shear_multiplier = 1
cross_multiplier = 1 * shear_multiplier

gFs = []
gFs.append(sa * wa * hf)
gFs.append(sb * wb * hf)
gF_cross_von_mises_a = cross_multiplier * sa * 2 * va * hc / np.sqrt(wb * wb / va * va + 3)
gF_cross_von_mises_b = cross_multiplier * sb * 2 * vb * hc / np.sqrt(wa * wa / vb * vb + 3)
cross_gFs = [gF_cross_von_mises_a, gF_cross_von_mises_b]
gFs.append(np.maximum(gF_cross_von_mises_a, gF_cross_von_mises_b))
# gFs.append(gF_cross_von_mises_a)
# gFs.append(gF_cross_von_mises_b)
gFs.append(taz * 2 * va * wa * shear_multiplier)
gFs.append(tbz * 2 * vb * wb * shear_multiplier)

minF = np.minimum.reduce(gFs)
stress = minF / ((wa + wb) * (hf + hc))
F = stress * (wa + wb) * (hf + hc)

gs = []
gs.append(1 - wa / 2 / wa_min)
gs.append(1 - wb / 2 / wb_min)
gs.append(1 - va / wa_min)
gs.append(1 - vb / wb_min)
gs.append(1 - hf / h_min)
gs.append(1 - hc / h_min)
gs.append((va + vb) / l_max - 1)
gs.append(1 - wa * hf * sa / F)  # tensile
gs.append(1 - wb * hf * sb / F)
cross_gs = [
    1 - 2 * hc / F * va * sa / np.sqrt(wb * wb / va / va + 3) * cross_multiplier,
    1 - 2 * hc / F * vb * sb / np.sqrt(wa * wa / vb / vb + 3) * cross_multiplier]
gs.append(np.minimum(cross_gs[0], cross_gs[1]))  # cross
gs.append(1 - 2 * va * wa * taz / F * shear_multiplier)  # z shear
gs.append(1 - 2 * vb * wb * tbz / F * shear_multiplier)

names = ['wa', 'wb', 'va', 'vb', 'hf', 'hc', 'design', 'tension a', 'tension b', 'cross', 'shear Z a', 'shear Z b']
cross_gs_names = ['shear/bend a', 'shear/bend b']

F = stress * (wa + wb) * (hf + hc)

best_idx = np.unravel_index(np.argmax(stress[:, :, :]), stress[:, :, :].shape)

print(f"best: stress: {stress[best_idx]}")
print(f" wb={wb[best_idx]:.2f}; va={va[best_idx]:.2f}; lmax={l_max:.2f}; hf={hf[best_idx]:.2f}; ")
print(f" wa={wa[best_idx]:.2f}; vb={vb[best_idx]:.2f}; hc={hc[best_idx]:.2f}; F={F[best_idx]:.2f}")

constraints = np.maximum.reduce(gs)
for i, g in enumerate(gs):
    if g.max() > .001:
        print(f"constraint g{i} is violated! : {g.max():.4f}")


print("\n== Active constraints: ", end="")
for i, g in enumerate(gs):
    if abs(g[best_idx]) < .01:
        print(names[i], end=", ")
print("")

print("-- cross beam constraints: ", end="")
for i, g in enumerate(cross_gs):
    if abs(g[best_idx]) < .01:
        print(cross_gs_names[i], end=", ")
print("\n")

# fig, axs = plt.subplots(2, 2, figsize=(10, 6))
# axs[0, 0].scatter(l.reshape(-1), FEM_stress.reshape(-1), c=colors)
# axs[0, 0].set_title('Total length')
# plt.show()



idx = best_idx # (2, 2, 2, 2)

fig, ax = plt.subplots(1, 1, figsize=(10, 6), subplot_kw={'projection': '3d'})
ax.plot_surface(hf[:, idx[1], :], wb[:, idx[1], :], stress[:, idx[1], :], edgecolor='none')
# ax.plot_surface(hf[idx[0], :, :], wb[idx[0], :, :], stress[idx[0], :, :], edgecolor='none')

plt.show()