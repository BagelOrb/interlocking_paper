import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib._color_data as mcd

from math import sqrt


bending = 0.0
cross_shear_force_ratio_shift = 1
combine_tensile_and_z_shear = False
full_stress_on_z = False
combine_z_shear_and_cross_shear = True

compare_to_FEM = False
broken_optimum = False

softmin_all_constraints = False
silicone = False

use_z_tensile_data = not compare_to_FEM
z_shear_stress_cross_beam_inclusion_a = float(full_stress_on_z or broken_optimum)
z_shear_stress_cross_beam_inclusion_b = float(full_stress_on_z or broken_optimum)
cross_beam_force_pillar_inclusion_a = float(broken_optimum)
cross_beam_force_pillar_inclusion_b = float(broken_optimum)
apply_cross_a = not broken_optimum
apply_cross_b = broken_optimum
use_z_shear_a = True  # not broken_optimum
use_z_shear_b = True

show_results = True
plot_legend = False

hf_sampling_multiplier = 10

lmax_val = 3.6
zoom_on_optimum = False
zoom = .05  # mm
optimum = (1.49, 0.35 / lmax_val, 0.5)

strain_a = 3.5 / 100
strain_b = 29 / 100

sa = 47
sb = 10.5
saz = 33
sbz = 10.6

if silicone:
    sb = sbz = 2.5

if not use_z_tensile_data:
    saz = sa
    sbz = sb

ta = sa / sqrt(3)  # TODO: reformulate constraints in terms of new gFs variables and remove these variables here
tb = sb / sqrt(3)
taz = saz / sqrt(3)
tbz = sbz / sqrt(3)

line_w = .3
layer_thickness = .1


#


wa_min = line_w
wb_min = line_w
w_max = 8 * 2 * wa_min
l_max = 6 * 2 * wa_min

h_min = 2 * layer_thickness
h_max = 10 * h_min


Nhf = 13 * hf_sampling_multiplier # 13 is whole layer heights
Nlmax = 1
Nwb = Nva = 400

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
    if zoom_on_optimum:
        wbs = np.linspace(max(2 * wb_min, optimum[0] - zoom), optimum[0] + zoom, Nwb)
        rys = np.linspace(optimum[1] - zoom, optimum[1] + zoom, Nva)
    else:
        wbs = np.linspace(2 * wb_min, w_max - 2 * wa_min, Nwb)
        rys = np.linspace(0.01, .99, Nva)
    hfs = np.linspace(h_min, h_max - h_min, Nhf)
    lmaxs = np.linspace(lmax_val if Nlmax == 1 else 3.6, lmax_val if Nlmax == 1 else 1.8, Nlmax)
    if silicone:
        lmaxs = np.linspace(0.8, 0.4, Nlmax)
    hf, lmax, wb, ry = np.meshgrid(hfs, lmaxs, wbs, rys, indexing='ij')
    va = ry * lmax

    if Nhf < 30:
        print(f"hfs: {hfs}")

# FEM_stress[np.isnan(FEM_stress)] = 0
# print( np.any(np.isnan( FEM_stress)))

wa = np.full(shape, 2 * wa_min)
hc = np.full(shape, h_min * (1 + float(silicone)))
vb = lmax - va
l = va + vb
w = wa + wb

def sq(x):
    return x * x

dya = (wb + wa * cross_beam_force_pillar_inclusion_a) / w
dyb = (wa + wb * cross_beam_force_pillar_inclusion_b) / w
dza = (wa + wb * z_shear_stress_cross_beam_inclusion_a) / w
dzb = (wb + wa * z_shear_stress_cross_beam_inclusion_b) / w

s11_a = 1 / (wa * hf)  # tensile
s11_b = 1 / (wb * hf)
s31_a = dya / (2 * va * hc)  # cross shear
s31_b = dyb / (2 * vb * hc)
s22_a = bending * s31_a * wb / va  # bending
s22_b = bending * s31_b * wa / vb
s12_a = dza / (2 * va * wa) * sa / saz  # Z shear
s12_b = dzb / (2 * vb * wb) * sb / sbz

if cross_shear_force_ratio_shift != 1:
    new_ratio = np.power(wa / w, cross_shear_force_ratio_shift)
    s31_a = s31_a / (wa / w) * new_ratio
    s31_b = s31_b / (wb / w) * (1 - new_ratio)

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
if use_z_shear_a and (not combine_z_shear_and_cross_shear or not apply_cross_a):
    gFs['Z shear a'] = sa / s12_a / sqrt(3)
if use_z_shear_b and (not combine_z_shear_and_cross_shear or not apply_cross_b):
    gFs['Z shear b'] = sb / s12_b / sqrt(3)
if combine_tensile_and_z_shear:
    gFs['tensile and z shear a'] = sa / np.sqrt(combined_von_mises_tensile_z_shear_a)
    gFs['tensile and z shear b'] = sb / np.sqrt(combined_von_mises_tensile_z_shear_b)
if combine_z_shear_and_cross_shear:
    z_shear_and_cross_a = sa / np.sqrt(combined_von_mises_z_shear_cross_a)
    z_shear_and_cross_b = sb / np.sqrt(combined_von_mises_z_shear_cross_b)
    if apply_cross_a and apply_cross_b and use_z_shear_b:
        gFs['z shear and cross a+b'] = np.maximum(z_shear_and_cross_a, z_shear_and_cross_b)
    else:
        if apply_cross_a and use_z_shear_a:
            gFs['z shear and cross a'] = z_shear_and_cross_a
        if apply_cross_b and use_z_shear_b:
            gFs['z shear and cross b'] = z_shear_and_cross_b
if apply_cross_a and apply_cross_b and not combine_z_shear_and_cross_shear:
    if bending == 0:
        cross_shear_a = sa / s31_a / sqrt(3)
        cross_shear_b = sb / s31_b / sqrt(3)
        gFs['cross shear a+b'] = np.maximum(cross_shear_a, cross_shear_b)
    else:
        cross_shear_and_cross_bending_a = sa / np.sqrt(combined_von_mises_cross_a)
        cross_shear_and_cross_bending_b = sa / np.sqrt(combined_von_mises_cross_b)
        gFs['cross shear and cross bending a+b'] = np.maximum(cross_shear_and_cross_bending_a, cross_shear_and_cross_bending_b)
elif apply_cross_a and (not combine_z_shear_and_cross_shear or not use_z_shear_b):
    gFs['cross shear a'] = sa / s31_a / sqrt(3)
    if bending > 0:
        gFs['cross shear and cross bending a'] = sa / np.sqrt(combined_von_mises_cross_a)
elif apply_cross_b and (not combine_z_shear_and_cross_shear or not use_z_shear_a):
    gFs['cross shear b'] = sb / s31_b / sqrt(3)
    if bending > 0:
        gFs['cross shear and cross bending b'] = sb / np.sqrt(combined_von_mises_cross_b)

gFs_mixed = np.zeros(shape)
if softmin_all_constraints:
    P = -10
    gFs_powered = []
    for gF in gFs.values():
        gFs_powered.append(np.power(gF, np.full(shape, P)))
    gFs_mixed = np.power(np.sum(gFs_powered, axis=0), np.full(shape, 1/P))


gMs = {}
gMs['wa'] = 1 - wa / 2 / wa_min
gMs['wb'] = 1 - wb / 2 / wb_min
gMs['va'] = 1 - va / wa_min
gMs['vb'] = 1 - vb / wb_min
gMs['hf'] = 1 - hf / h_min
gMs['hc'] = 1 - hc / h_min


print("Used constraints: ", end="")
for name, gF in gFs.items():
    print(name + ", ", end="")
print("\n")

minF = gFs_mixed if softmin_all_constraints else np.minimum.reduce(list(gFs.values()))
stress = minF / ((wa + wb) * (hf + hc))

valid_manufacturing_constraints_multiplier = ((np.maximum.reduce(list(gMs.values())) <= 0.0001) + .00001) / 1 + .00001
stress = stress * valid_manufacturing_constraints_multiplier
FEM_stress = FEM_stress * valid_manufacturing_constraints_multiplier

F = stress * (wa + wb) * (hf + hc)
if compare_to_FEM:
    F_FEM = FEM_stress * (wa + wb) * (hf + hc)

# F * sqrt(term) < s
# 1 - s / (F * sqrt(term)) < 0
cross_gs = [
    1 - sa / (F * np.sqrt(combined_von_mises_cross_a)),
    1 - sb / (F * np.sqrt(combined_von_mises_cross_b))]
cross_gs_names = ['shear/bend a', 'shear/bend b']

gs = {}
gs['design'] = (va + vb) / l_max - 1
gs['tensile a'] = 1 - wa * hf * sa / F
gs['tensile b'] = 1 - wb * hf * sb / F
if use_z_shear_a:
    gs['z shear a'] = 1 - 2 * va * wa * taz / (F * dza)
if use_z_shear_b:
    gs['z shear b'] = 1 - 2 * vb * wb * tbz / (F * dzb)
if combine_tensile_and_z_shear:
    gs['tensile and z shear a'] = F * F * combined_von_mises_tensile_z_shear_a / sq(sa) - 1
    if use_z_shear_b:
        gs['tensile and z shear b'] = F * F * combined_von_mises_tensile_z_shear_b / sq(sb) - 1
if combine_z_shear_and_cross_shear:
    z_shear_and_cross_a = F * F * combined_von_mises_z_shear_cross_a / sq(sa) - 1
    z_shear_and_cross_b = F * F * combined_von_mises_z_shear_cross_b / sq(sb) - 1
    if apply_cross_a and apply_cross_b and use_z_shear_b:
        gs['z shear and cross a+b'] = np.minimum(z_shear_and_cross_a, z_shear_and_cross_b)
    else:
        gs['z shear and cross a'] = z_shear_and_cross_a
        if use_z_shear_b:
            gs['z shear and cross b'] = z_shear_and_cross_b

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
        if abs(g[best_idx]) < .01:
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
            print(f"{name}:  {np.average(ratios):.1%}, stdev: {np.std(ratios):.1%}")
        else:
            print(f"{name}")

    print("- cross beam sub failure modes -")
    for i, gF in enumerate(cross_gFs):
        ratios = prediction_ratio[minF == gF]
        if ratios.size > 0:
            print(f"g{i} {cross_gs_names[i]}:  {np.average(ratios):.1%}, stdev: {np.std(ratios):.1%}")
        else:
            print(f"g{i} {cross_gs_names[i]}")

constraints = np.maximum.reduce(gs)
for name, g in gs.items():
    best_idx = np.unravel_index(np.argmax(stress[:, :, :, :]), shape)
    if g.max() > .001:
        print(f"!!! > {name} constraint is violated! : {g.max():.4f}")
    else:
        print(f"{name}: {g[best_idx]}")

#


#

if show_results:
    colormap = {'tensile a': "green",
                'tensile b': "red",
                'Z shear a': "lime",
                'Z shear b': "tomato",
                'tensile and z shear a': "green",
                'tensile and z shear b': "red",
                'z shear and cross a+b': "magenta",
                'z shear and cross a': "darkgreen",
                'z shear and cross b': "maroon",
                'cross shear a+b': "blue",
                'cross shear and cross bending a+b': "aquamarine",
                'cross shear a': "yellowgreen",
                'cross shear b': "orange",
                'cross shear and cross bending a': "chartreuse",
                'cross shear and cross bending b': "pink"}

    failure_mode_colors = np.full(shape, "black", dtype=object)
    for name, gF in gFs.items():
        failure_mode_colors[minF == gF] = mcd.XKCD_COLORS["xkcd:"+colormap[name]] + "c0"


    #


    def plotTwo(ax, X, Y, Z1, Z2, col1=None, col2=None):
        if col1 is None or softmin_all_constraints:
            col1 = np.full(Z1.shape, "#00ff00c0", dtype=object)
        if col2 is None:
            col2 = np.full(Z2.shape, "#999999c0", dtype=object)
        if compare_to_FEM:
            col1[-1, :] = "#ffffff00"
            ax.plot_surface(np.append(X, X, axis=0),
                            np.append(Y, Y, axis=0),
                            np.append(Z1, Z2, axis=0),
                            facecolors=np.append(col1, col2, axis=0),
                            edgecolor='none')
        else:
            ax.plot_surface(X, Y, Z1, facecolors=col1, edgecolor='none')

    if not compare_to_FEM:
        FEM_stress = stress

    best_FEM_idx = np.unravel_index(np.argmax(FEM_stress), FEM_stress.shape)
    idx = best_FEM_idx

    # hf, lmax, wb, va
    l = 0

    fig, ax = plt.subplots(1, 3, subplot_kw={'projection': '3d'})
    plotTwo(ax[0], wb[idx[0], l, :, :], va[idx[0], l, :, :], stress[idx[0], l, :, :], FEM_stress[idx[0], l, :, :], failure_mode_colors[idx[0], l, :, :])
    ax[0].set(xlabel='wb', ylabel='va')
    plotTwo(ax[1], hf[:, l, idx[2], :], va[:, l, idx[2], :], stress[:, l, idx[2], :], FEM_stress[:, l, idx[2], :], failure_mode_colors[:, l, idx[2], :])
    ax[1].set(xlabel='hf', ylabel='va')
    plotTwo(ax[2], hf[:, l, :, idx[3]], wb[:, l, :, idx[3]], stress[:, l, :, idx[3]], FEM_stress[:, l, :, idx[3]], failure_mode_colors[:, l, :, idx[3]])
    ax[2].set(xlabel='hf', ylabel='wb')

    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')

    if plot_legend:
        fig, ax = plt.subplots()
        legend_elements = []
        for name, color in colormap.items():
            legend_elements.append(Patch(facecolor=color, edgecolor='none', label=name))
        ax.legend(handles=legend_elements)

    # plt.savefig("C:\\Users\\t.kuipers\\OneDrive - Ultimaker B.V\\Documents\\PhD\\interlocking_project\\paper\\paper_git\\analytical_optimization\\ana_correct_z_yield_ORcross.svg")

    plt.show()

