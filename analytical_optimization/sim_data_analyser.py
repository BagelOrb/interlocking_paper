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
h_max = 6 * h_min

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

# FEM_stress[np.isnan(FEM_stress)] = 0
# print( np.any(np.isnan( FEM_stress)))

wa = np.full(shape, 2 * wa_min)
hc = np.full(shape, h_min)
vb = lmax - va
l = va + vb

F_FEM = FEM_stress * (wa + wb) * (hf + hc)

w = wa + wb

shear_multiplier = 1
cross_multiplier = 1 * shear_multiplier
bending = 0
z_shear_stress_cross_beam_inclusion = 1
apply_cross_b = False
tensile_in_von_mises = 0
use_combined_von_mises = False


def sq(x):
    return x * x


s12_a = 1 / (2 * va * w)  # Z shear
s31_a = s12_a * wb / hc  # cross shear
s22_a = bending * s31_a * wb / va  # bending
s11_a = tensile_in_von_mises * 1 / (wa * hf)  # tensile

s12_b = 1 / (2 * vb * w)
s31_b = s12_b * wa / hc
s11_b = tensile_in_von_mises * 1 / (wb * hf)
s22_b = bending * s31_b * wa / vb

term_a = sq(s11_a) +                        sq(s22_a) - s11_a * s22_a + 3 * (                       sq(s31_a) + sq(s12_a))
term_b = sq(s11_b) + float(apply_cross_b) * sq(s22_b) - s11_b * s22_b + 3 * (float(apply_cross_b) * sq(s31_b) + sq(s12_b))

gF_cross_von_mises_a = 1 / (wb / w) * cross_multiplier * sa * 2 * va * hc / np.sqrt(bending * wb * wb / va * va + 3)
gF_cross_von_mises_b = 1 / (wa / w) * cross_multiplier * sb * 2 * vb * hc / np.sqrt(bending * wa * wa / vb * vb + 3)
cross_gFs = [gF_cross_von_mises_a, gF_cross_von_mises_b]

gFs = []
if not use_combined_von_mises or tensile_in_von_mises == 0:
    gFs.append(sa * wa * hf)
    gFs.append(sb * wb * hf)
if use_combined_von_mises:
    gFs.append(np.sqrt(sa * sa / term_a))
    gFs.append(np.sqrt(sb * sa / term_b))
else:
    if apply_cross_b:
        gFs.append(np.maximum(gF_cross_von_mises_a, gF_cross_von_mises_b))
    else:
        gFs.append(gF_cross_von_mises_a)
    gFs.append(1 / ((wa + z_shear_stress_cross_beam_inclusion * wb) / w) * taz * 2 * va * wa * shear_multiplier)
    gFs.append(1 / ((wb + z_shear_stress_cross_beam_inclusion * wa) / w) * tbz * 2 * vb * wb * shear_multiplier)

minF = np.minimum.reduce(gFs)
stress = minF / ((wa + wb) * (hf + hc))
F = stress * (wa + wb) * (hf + hc)

cross_gs = [
    1 - 2 * hc / (F * wb / (wa + wb)) * va * sa / np.sqrt(bending * wb * wb / va / va + 3) * cross_multiplier,
    1 - 2 * hc / (F * wa / (wa + wb)) * vb * sb / np.sqrt(bending * wa * wa / vb / vb + 3) * cross_multiplier]

gs = []
gs.append(1 - wa / 2 / wa_min)
gs.append(1 - wb / 2 / wb_min)
gs.append(1 - va / wa_min)
gs.append(1 - vb / wb_min)
gs.append(1 - hf / h_min)
gs.append(1 - hc / h_min)
gs.append((va + vb) / l_max - 1)
if not use_combined_von_mises or tensile_in_von_mises == 0:
    gs.append(1 - wa * hf * sa / F)  # tensile
    gs.append(1 - wb * hf * sb / F)
if use_combined_von_mises:
    gs.append(F * F * term_a / sq(sa) - 1)
    gs.append(F * F * term_b / sq(sb) - 1)
else:
    if apply_cross_b:
        gs.append(np.minimum(cross_gs[0], cross_gs[1]))  # cross
    else:
        gs.append(cross_gs[0])
    gs.append(1 - 2 * va * wa * taz / (F * (wa + z_shear_stress_cross_beam_inclusion * wb) / w) * shear_multiplier)  # z shear
    gs.append(1 - 2 * vb * wb * tbz / (F * (wa + z_shear_stress_cross_beam_inclusion * wb) / w) * shear_multiplier)


names = ['wa', 'wb', 'va', 'vb', 'hf', 'hc', 'design', 'tension a', 'tension b', 'cross', 'shear Z a', 'shear Z b', 'combined a', 'combined b']
cross_gs_names = ['shear/bend a', 'shear/bend b']

for l in range(Nlmax):
    print(f"\n\n\n===  L_max: {lmax[0, l, 0, 0]}\n\n")

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

    print("\n== prediction ratio per failure mode ==")
    prediction_ratio = stress / FEM_stress
    for i, gF in enumerate(gFs):
        ratios = prediction_ratio[:, l, :, :][minF[:, l, :, :] == gF[:, l, :, :]]
        if ratios.size > 0:
            print(f"g{i} {names[i + 7]}:  {np.average(ratios):.3f}, stdev: {np.std(ratios):.3f}")
        else:
            print(f"g{i} {names[i + 7]}")

    print("- cross beam sub failure modes -")
    for i, gF in enumerate(cross_gFs):
        ratios = prediction_ratio[:, l, :, :][minF[:, l, :, :] == gF[:, l, :, :]]
        if ratios.size > 0:
            print(f"g{i} {cross_gs_names[i]}:  {np.average(ratios):.3f}, stdev: {np.std(ratios):.3f}")
        else:
            print(f"g{i} {cross_gs_names[i]}")

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
    col1 = np.full(Z1.shape, 'r', dtype='U50')
    col2 = np.full(Z2.shape, 'b', dtype='U50')
    ax.plot_surface(np.append(X, X, axis=0),
                    np.append(Y, Y, axis=0),
                    np.append(Z1, Z2, axis=0),
                    facecolors= np.append(col1, col2, axis=0),
                    edgecolor='none', alpha=.75)


idx = best_FEM_idx  # (2, 2, 2, 2)
# hf, lmax, wb, va

fig, ax = plt.subplots(1, 3, figsize=(10, 6), subplot_kw={'projection': '3d'})
plotTwo(ax[0], wb[idx[0], 0, :, :], va[idx[0], 0, :, :], stress[idx[0], 0, :, :], FEM_stress[idx[0], 0, :, :])
ax[0].set(xlabel='wb', ylabel='va')
plotTwo(ax[1], hf[:, 0, idx[2], :], va[:, 0, idx[2], :], stress[:, 0, idx[2], :], FEM_stress[:, 0, idx[2], :])
ax[1].set(xlabel='hf', ylabel='va')
plotTwo(ax[2], hf[:, 0, :, idx[3]], wb[:, 0, :, idx[3]], stress[:, 0, :, idx[3]], FEM_stress[:, 0, :, idx[3]])
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

mng = plt.get_current_fig_manager()
mng.full_screen_toggle()

plt.show()

# print(np.maximum.reduce(gs))
