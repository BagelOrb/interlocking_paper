import numpy as np

import matplotlib.pyplot as plt
from matplotlib import pyplot
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

from math import sqrt

strain_a = 3.5 / 100
strain_b = 29 / 100

sa = 47  # / (1 + strain_a)
sb = 10.5  # / (1 + strain_b)
saz = sa  # 33
sbz = sb  # 9.0

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
h_max = 10 * h_min

compare_to_FEM = False

Nhf = 90
Nlmax = 1
Nwb = Nva = 45

if compare_to_FEM:
    # reading file
    data_file = np.genfromtxt('IS-data-0527.csv', delimiter=',')
    data_file[0, 0] = 0.2  # why doesn't it read that cell?!?!

    Nhf = 5
    Nlmax = 4
    Nwb = Nva = 9

shape = (Nhf, Nlmax, Nwb, Nva)

FEM_stress = np.ndarray(shape=shape)
hf = np.ndarray(shape=shape)
lmax = np.ndarray(shape=shape)
wb = np.ndarray(shape=shape)
va = np.ndarray(shape=shape)

if compare_to_FEM:
    for hf_ in range(0, Nhf):
        hf[hf_, :, :, :] = data_file[(Nwb + 1) * hf_, 0]

    for lmax_ in range(0, Nlmax):
        lmax[:, lmax_, :, :] = data_file[0, 1 + (Nva + 1) * lmax_]

    for wb_ in range(0, Nwb):
        wb[:, :, wb_, :] = data_file[1 + wb_, 1]

    for va_ in range(0, Nva):
        for lmax_ in range(0, Nlmax):
            va[:, lmax_:, :, va_] = data_file[0, 2 + lmax_ * (Nva + 1) + va_]

    for hf_ in range(0, Nhf):
        for lmax_ in range(0, Nlmax):
            FEM_stress[hf_, lmax_, :, :] = data_file[1 + (Nwb + 1) * hf_:1 + (Nwb + 1) * hf_ + Nwb,
                                           2 + (Nva + 1) * lmax_:2 + (Nva + 1) * lmax_ + Nva]
else:
    wbs = np.linspace(2 * wb_min, w_max - 2 * wa_min, Nwb)
    rys = np.linspace(0.1, .9, Nva)
    hfs = np.linspace(h_min, h_max - h_min, Nhf)
    lmaxs = np.linspace(3.6, 1.8, Nlmax)
    hf, lmax, wb, ry = np.meshgrid(hfs, lmaxs, wbs, rys, indexing='ij')
    va = ry * lmax

    if Nhf < 30:
        print(f"hfs: {hfs}")

# FEM_stress[np.isnan(FEM_stress)] = 0
# print( np.any(np.isnan( FEM_stress)))

wa = np.full(shape, 2 * wa_min)
hc = np.full(shape, h_min)
vb = lmax - va
l = va + vb
w = wa + wb

if compare_to_FEM:
    F_FEM = FEM_stress * (wa + wb) * (hf + hc)

bending = 0.0
z_shear_stress_cross_beam_inclusion = 1
apply_cross_a = True
apply_cross_b = False  # TODO: fix applying OR for breakage of both beams
combine_tensile_and_z_shear = False
combine_z_shear_and_cross_shear = True


def sq(x):
    return x * x


assert (apply_cross_a)

s12_a = 1 / (2 * va * w)  # Z shear
s12_b = 1 / (2 * vb * w)
s31_a = s12_a * wb / hc  # cross shear
s31_b = s12_b * wa / hc
s22_a = bending * s31_a * wb / va  # bending
s22_b = bending * s31_b * wa / vb
s11_a = 1 / (wa * hf)  # tensile
s11_b = 1 / (wb * hf)

if z_shear_stress_cross_beam_inclusion > 0:
    s12_a = s12_a * (wa + z_shear_stress_cross_beam_inclusion * wb) / wa
    s12_b = s12_b * (wb + z_shear_stress_cross_beam_inclusion * wa) / wb

combined_von_mises_z_shear_cross_a = float(apply_cross_a) * sq(s22_a) + 3 * (float(apply_cross_a) * sq(s31_a) + sq(s12_a))
combined_von_mises_z_shear_cross_b = float(apply_cross_b) * sq(s22_b) + 3 * (float(apply_cross_b) * sq(s31_b) + sq(s12_b))
combined_von_mises_tensile_z_shear_a = sq(s11_a) + 3 * sq(s12_a)
combined_von_mises_tensile_z_shear_b = sq(s11_b) + 3 * sq(s12_b)

combined_von_mises_cross_a = sq(s22_a) + 3 * sq(s31_a)
combined_von_mises_cross_b = sq(s22_b) + 3 * sq(s31_b)

gF_cross_von_mises_a = 1 / (wb / w) * sa * 2 * va * hc / np.sqrt(sq(bending * wb / va) + 3)
gF_cross_von_mises_b = 1 / (wa / w) * sb * 2 * vb * hc / np.sqrt(sq(bending * wa / vb) + 3)
cross_gFs = [gF_cross_von_mises_a, gF_cross_von_mises_b]

gFs = {}
gFs['tensile a'] = sa / s11_a
gFs['tensile b'] = sb / s11_b
gFs['Z shear a'] = sa / s12_a / sqrt(3)
gFs['Z shear b'] = sb / s12_b / sqrt(3)
if combine_tensile_and_z_shear:
    gFs['tensile and z shear a'] = sa / np.sqrt(combined_von_mises_tensile_z_shear_a)
    gFs['tensile and z shear b'] = sb / np.sqrt(combined_von_mises_tensile_z_shear_b)
if combine_z_shear_and_cross_shear:
    gFs['z shear and cross a'] = sa / np.sqrt(combined_von_mises_z_shear_cross_a)
    gFs['z shear and cross b'] = sb / np.sqrt(combined_von_mises_z_shear_cross_b)
if apply_cross_a and apply_cross_b:
    if bending == 0:
        cross_shear_a = sa / s31_a / sqrt(3)
        cross_shear_b = sb / s31_b / sqrt(3)
        gFs['cross shear a+b'] = np.maximum(cross_shear_a, cross_shear_b)
    else:
        cross_shear_and_cross_bending_a = sa / np.sqrt(combined_von_mises_cross_a)
        cross_shear_and_cross_bending_b = sa / np.sqrt(combined_von_mises_cross_b)
        gFs['cross shear and cross bending a+b'] = np.maximum(cross_shear_and_cross_bending_a, cross_shear_and_cross_bending_b)
elif apply_cross_a:
    gFs['cross shear a'] = sa / s31_a / sqrt(3)
    if bending > 0:
        gFs['cross shear and cross bending a'] = sa / np.sqrt(combined_von_mises_cross_a)
elif apply_cross_b:
    gFs['cross shear b'] = sb / s31_b / sqrt(3)
    if bending > 0:
        gFs['cross shear and cross bending b'] = sb / np.sqrt(combined_von_mises_cross_b)


minF = np.minimum.reduce(list(gFs.values()))
stress = minF / ((wa + wb) * (hf + hc))
F = stress * (wa + wb) * (hf + hc)

cross_gs = [
    1 - 2 * hc / (F * wb / (wa + wb)) * va * sa / np.sqrt(bending * wb * wb / va / va + 3),
    1 - 2 * hc / (F * wa / (wa + wb)) * vb * sb / np.sqrt(bending * wa * wa / vb / vb + 3)]
cross_gs_names = ['shear/bend a', 'shear/bend b']

gs = {}
gs['wa'] = 1 - wa / 2 / wa_min
gs['wb'] = 1 - wb / 2 / wb_min
gs['va'] = 1 - va / wa_min
gs['vb'] = 1 - vb / wb_min
gs['hf'] = 1 - hf / h_min
gs['hc'] = 1 - hc / h_min
gs['design'] = (va + vb) / l_max - 1
gs['tensile a'] = 1 - wa * hf * sa / F
gs['tensile b'] = 1 - wb * hf * sb / F
gs['z shear a'] = 1 - 2 * va * wa * taz / (F * (wa + z_shear_stress_cross_beam_inclusion * wb) / w)
gs['z shear b'] = 1 - 2 * vb * wb * tbz / (F * (wa + z_shear_stress_cross_beam_inclusion * wb) / w)
if combine_tensile_and_z_shear:
    gs['tensile and z shear a'] = F * F * combined_von_mises_tensile_z_shear_a / sq(sa) - 1
    gs['tensile and z shear b'] = F * F * combined_von_mises_tensile_z_shear_b / sq(sb) - 1
if combine_z_shear_and_cross_shear:
    gs['z shear and cross a'] = F * F * combined_von_mises_z_shear_cross_a / sq(sa) - 1
    gs['z shear and cross b'] = F * F * combined_von_mises_z_shear_cross_b / sq(sb) - 1

if apply_cross_a and apply_cross_b:
    gs['cross a+b'] = np.minimum(cross_gs[0], cross_gs[1])  # cross
elif apply_cross_a:
    gs['cross a'] = cross_gs[0]
elif apply_cross_b:
    gs['cross b'] = cross_gs[1]

for l in range(Nlmax):
    print(f"\n\n\n===  L_max: {lmax[0, l, 0, 0]}\n\n")

    if compare_to_FEM:
        best_FEM_idx = np.unravel_index(np.argmax(FEM_stress[:, l, :, :]), FEM_stress[:, l, :, :].shape)
        best_FEM_idx = (best_FEM_idx[0], l, best_FEM_idx[1], best_FEM_idx[2])
        print(f"best FEM: stess: {FEM_stress[best_FEM_idx]}")
        print(f" wb={wb[best_FEM_idx]:.2f}; va={va[best_FEM_idx]:.2f}; lmax={lmax[best_FEM_idx]:.2f}; hf={hf[best_FEM_idx]:.2f}; ")
        print(f" wa={wa[best_FEM_idx]:.2f}; vb={vb[best_FEM_idx]:.2f}; hc={hc[best_FEM_idx]:.2f}; F={F[best_FEM_idx]:.2f}")
        print("")

    best_idx = np.unravel_index(np.argmax(stress[:, l, :, :]), stress[:, l, :, :].shape)
    best_idx = (best_idx[0], l, best_idx[1], best_idx[2])
    print(f"best ana: stress: {stress[best_idx]}")
    print(f" wb={wb[best_idx]:.2f}; va={va[best_idx]:.2f}; lmax={lmax[best_idx]:.2f}; hf={hf[best_idx]:.2f}; ")
    print(f" wa={wa[best_idx]:.2f}; vb={vb[best_idx]:.2f}; hc={hc[best_idx]:.2f}; F={F[best_idx]:.2f}")

    print("\n== Active constraints: ", end="")
    for name, g in gs.items():
        if abs(g[best_idx]) < .001:
            print(name, end=",")
    print("")

    print("-- cross beam constraints: ", end="")
    for i, g in enumerate(cross_gs):
        if abs(g[best_idx]) < .001:
            print(cross_gs_names[i], end=",")
    print("\n")

if compare_to_FEM:
    prediction_ratio = stress / FEM_stress
    print(f"\n== prediction ratio : {np.average(prediction_ratio):.1%}, stdev: {np.std(prediction_ratio):.1%}")
    print("== prediction ratio per failure mode ==")
    for name, gF in gFs.items():
        ratios = prediction_ratio[minF == gF]
        if ratios.size > 0:
            print(f"g{i} {name}:  {np.average(ratios):.1%}, stdev: {np.std(ratios):.1%}")
        else:
            print(f"g{i} {name}")

    print("- cross beam sub failure modes -")
    for i, gF in enumerate(cross_gFs):
        ratios = prediction_ratio[minF == gF]
        if ratios.size > 0:
            print(f"g{i} {cross_gs_names[i]}:  {np.average(ratios):.1%}, stdev: {np.std(ratios):.1%}")
        else:
            print(f"g{i} {cross_gs_names[i]}")

constraints = np.maximum.reduce(gs)
for name, g in gs.items():
    if g.max() > .001:
        print(f"{name} constraint is violated! : {g.max():.4f}")


#


#


colors = np.zeros((Nhf * Nlmax * Nwb * Nva, 3))
colors[:, 0] = wb.reshape(-1) / wb.max()
colors[:, 1] = va.reshape(-1) / va.max()
colors[:, 2] = hf.reshape(-1) / hf.max()
# colors[:, 2] = va.reshape(-1) / va.max()

if False:
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


def plotTwo(ax, X, Y, Z1, Z2):
    if compare_to_FEM:
        col1 = np.full(Z1.shape, "#00ff00d0", dtype=object)
        col2 = np.full(Z2.shape, "#ff0000d0", dtype=object)
        col1[-1, :] = "#ffffff00"
        ax.plot_surface(np.append(X, X, axis=0),
                        np.append(Y, Y, axis=0),
                        np.append(Z1, Z2, axis=0),
                        facecolors=np.append(col1, col2, axis=0),
                        edgecolor='none')
    else:
        ax.plot_surface(X, Y, Z1, color='g', edgecolor='none')



if not compare_to_FEM:
    FEM_stress = stress

best_FEM_idx = np.unravel_index(np.argmax(FEM_stress), FEM_stress.shape)
idx = best_FEM_idx

# hf, lmax, wb, va
l = 0

fig, ax = plt.subplots(1, 3, figsize=(10, 6), subplot_kw={'projection': '3d'})
plotTwo(ax[0], wb[idx[0], l, :, :], va[idx[0], l, :, :], stress[idx[0], l, :, :], FEM_stress[idx[0], l, :, :])
ax[0].set(xlabel='wb', ylabel='va')
plotTwo(ax[1], hf[:, l, idx[2], :], va[:, l, idx[2], :], stress[:, l, idx[2], :], FEM_stress[:, l, idx[2], :])
ax[1].set(xlabel='hf', ylabel='va')
plotTwo(ax[2], hf[:, l, :, idx[3]], wb[:, l, :, idx[3]], stress[:, l, :, idx[3]], FEM_stress[:, l, :, idx[3]])
ax[2].set(xlabel='hf', ylabel='wb')

if False:
    fig, axs = plt.subplots(1, 3, figsize=(10, 6), subplot_kw={'projection': '3d'})
    axs[0].plot_surface(hf[:, idx[1], :, idx[3]], wb[:, idx[1], :, idx[3]], FEM_stress[:, idx[1], :, idx[3]], edgecolor='none')
    axs[0].set_title('FEM')
    axs[1].plot_surface(hf[:, idx[1], :, idx[3]], wb[:, idx[1], :, idx[3]], FEM_stress[:, idx[1], :, idx[3]], edgecolor='none', alpha=.75)
    axs[1].plot_surface(hf[:, idx[1], :, idx[3]], wb[:, idx[1], :, idx[3]], stress[:, idx[1], :, idx[3]], edgecolor='none', alpha=.75)
    axs[1].set_title('combined')
    axs[2].plot_surface(hf[:, idx[1], :, idx[3]], wb[:, idx[1], :, idx[3]], stress[:, idx[1], :, idx[3]], edgecolor='none')
    axs[2].set_title('analytical')

wm = plt.get_current_fig_manager()
wm.window.state('zoomed')

plt.show()

# print(np.maximum.reduce(gs))
