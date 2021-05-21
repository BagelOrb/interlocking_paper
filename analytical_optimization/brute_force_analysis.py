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
gFs.append(.5 * np.maximum(sa * va, sb * vb) * 4 / 3 * hc)
# gFs.append(.5 * sb * 4 / 3 * vb * hc)
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
    f"rz={rz[best_idx]:.4g}; rw={rw[best_idx]:.4g}; rx={rx[best_idx]:.4g}; ry={ry[best_idx]:.4g}; FF={F_analytic[best_idx]:.4g};\n"
    f"wa={wa[best_idx]:.4g}; wb={wb[best_idx]:.4g}; va={va[best_idx]:.4g}; vb={vb[best_idx]:.4g}; hf={hf[best_idx]:.4g}; hc={hc[best_idx]:.4g};")



# fig, axs = plt.subplots(2, 2, figsize=(10, 6))
# axs[0, 0].scatter(l.reshape(-1), FEM_stress.reshape(-1), c=colors)
# axs[0, 0].set_title('Total length')
# plt.show()
