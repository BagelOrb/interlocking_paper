import numpy as np

import matplotlib.pyplot as plt
from math import sqrt

sa = 47
sb = 10.5
saz = 33
sbz = 9.0

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
    wb_opt = 1.8848
    va_opt = 1.9667
    hf_opt = 0.6202
    d = .1
    wbs = np.linspace(wb_opt - d, wb_opt + d, Nwb)
    vas = np.linspace(va_opt - d, va_opt + d, Nva)
    hfs = np.linspace(hf_opt - d, hf_opt + d, Nhf)

wb, va, hf = np.meshgrid(wbs, vas, hfs, indexing='ij')

l = np.full(shape, l_max)
wa = np.full(shape, 2 * wa_min)
hc = np.full(shape, h_min)
vb = l - va


gFs = []
gFs.append(sa * wa * hf)
gFs.append(sb * wb * hf)
gF_shear_a = 4 / 3 * hc / sqrt(3) * va * sa
gF_shear_b = 4 / 3 * hc / sqrt(3) * vb * sb
gF_bend_a = 2 * hc / wb * va * va * sa
gF_bend_b = 2 * hc / wa * vb * vb * sb
gF_cross_a = np.minimum(gF_shear_a, gF_bend_a)
gF_cross_b = np.minimum(gF_shear_b, gF_bend_b)
gFs.append(np.maximum(gF_cross_a, gF_cross_b))
# gFs.append(gF_cross_a)
gFs.append(taz * 4 / 3 * va * wa)
gFs.append(tbz * 4 / 3 * vb * wb)

minF = np.minimum.reduce(gFs)
stress = minF / ((wa + wb) * (hf + hc))
best_idx = np.unravel_index(np.argmax(stress), stress.shape)

F = stress * (wa + wb) * (hf + hc)
print(f"best: stress: {stress[best_idx]}")
print(f" wb={wb[best_idx]:.2f}; va={va[best_idx]:.2f}; lmax={l_max:.2f}; hf={hf[best_idx]:.2f}; ")
print(f" wa={wa[best_idx]:.2f}; vb={vb[best_idx]:.2f}; hc={hc[best_idx]:.2f}; F={F[best_idx]:.2f}")


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


names = ['wa', 'wb', 'va', 'vb', 'hf', 'hc', 'design', 'tension a', 'tension b', 'cross', 'shear Z a', 'shear Z b']
cross_gs_names = ['shear a', 'shear b', 'bend a', 'bend b']

constraints = np.maximum.reduce(gs)
assert (constraints.max() <= 0.00000001)

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
