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

rz = hf / hc
rw = va / wa
rx = wb / wa
ry = vb / va

gFs = []
gFs.append(sa * wa * hf)
gFs.append(sb * wb * hf)
gFs.append(4 / 3 * hc / np.minimum(np.maximum(2, wb / va) / (va * sa), np.maximum(2, wa / vb) / (vb * sb)))
gFs.append(taz * 4 / 3 * va * wa)
gFs.append(tbz * 4 / 3 * vb * wb)
gFs.append(sa * 4 / 3 * va * va * hc / wb)
gFs.append(sb * 4 / 3 * vb * vb * hc / wa)

minF = np.minimum.reduce(gFs)
stress = minF / ((wa + wb) * (hf + hc))
best_idx = np.unravel_index(np.argmax(stress), stress.shape)

F = stress * (wa + wb) * (hf + hc)
print(
    f"best: stress: {stress[best_idx]}\n"
    f"rz={rz[best_idx]:.4f}; rw={rw[best_idx]:.4f}; rx={rx[best_idx]:.4f}; ry={ry[best_idx]:.4f}; FF={F[best_idx]:.4f};\n"
    f"wa={wa[best_idx]:.4f}; wb={wb[best_idx]:.4f}; va={va[best_idx]:.4f}; vb={vb[best_idx]:.4f}; hf={hf[best_idx]:.4f}; hc={hc[best_idx]:.4f};")

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
gs.append(3 * F / (4 * vb * wb * tbz) - 1)

cross_gs = [
    3 * F / (4 * va * hc * ta) - 1  # shear
    , 3 * F / (4 * vb * hc * tb) - 1
    , 3 * F * wb / (4 * va * va * hc * sa) - 1  # bending
    , 3 * F * wa / (4 * vb * vb * hc * sb) - 1]

names = ['wa', 'wb', 'va', 'vb', 'hf', 'hc', 'design', 'tension a', 'tension b', 'cross', 'shear Z a', 'shear Z b']

constraints = np.maximum.reduce(gs)
assert (constraints.max() <= 0.00000001)

for i, g in enumerate(gs):
    print(names[i], f"{g[best_idx]:.4f}")

for g in cross_gs:
    print(f"{g[best_idx]:.4f}")
#

# fig, axs = plt.subplots(2, 2, figsize=(10, 6))
# axs[0, 0].scatter(l.reshape(-1), FEM_stress.reshape(-1), c=colors)
# axs[0, 0].set_title('Total length')
# plt.show()
