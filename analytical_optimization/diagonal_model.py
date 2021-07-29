import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib._color_data as mcd
import matplotlib.colors

from scipy import interpolate

from math import sqrt
from math import pi

big = True
only_one = False

compare_to_FEM = True
replace_by_FEM = False
make_comparison_plot = True

data_file = np.genfromtxt('diagonal_sim_results.csv', delimiter=';')
FEM_stress = data_file[1:, 1:].T


Ls = np.linspace(.6, 3.6, 200)
wbs = np.linspace(.6, 4.8, 200)
if not big:
    Ls = np.asarray([3.6])
    wbs = np.linspace(.6, 3.6, 6)
if compare_to_FEM:
    Ls = np.linspace(1.8, 3.6, 4)
    wbs = np.linspace(.6, 1.8, 5)
L, wb = np.meshgrid(Ls, wbs)
if only_one:
    L = np.asarray([3.6])
    wb = np.asarray([.6])
shape = L.shape

wmin = 0.3
r = .15
wa = 2 * wmin
h = .2

sa = 47
sb = 10.5
saz = 33
sbz = 10.6
if compare_to_FEM: saz = sa; sbz = sb

w = wa + wb


M = L - 2*r
div = np.sqrt(np.maximum(0.0001, 4*M**2 - w**2))
div[np.logical_not(np.isfinite(div))] = 999999
d = 2 * M * w / div
if only_one: print(f"{d=}")

Aza = .5 * pi * r**2 + r * (wa - 2*r)*d/w + d*M*(wa/w)**2
Azb = .5 * pi * r**2 + r * (wb - 2*r)*d/w + d*M*(wb/w)**2

div = np.sqrt((w/d)**2 + 3*(d/2/M)**2)
div[np.logical_not(np.isfinite(div))] = 99999
if only_one: print(f"{div=}")
gFs = {}
gFs["tensile a"] = 2 * sa * wa * h / div
gFs["tensile b"] = 2 * sb * wb * h / div
gFs["Z shear a"] = 2 * saz * Aza / sqrt(3)
gFs["Z shear b"] = 2 * sbz * Azb / sqrt(3)


minF = np.minimum.reduce(list(gFs.values()))
np.where(np.isnan(minF), minF, 0)

stress = minF / d / (2 * h)

if replace_by_FEM:
    stress = FEM_stress

if not big:
    print(f"{stress=}")


best_wb_idx = np.unravel_index(np.argmax(stress), stress.shape)
print(f"max stress: {stress.max()} at wb: {wb[best_wb_idx]}")

colormap = {"tensile a": "green",
            "tensile b": "red",
            "Z shear a": "orange",
            "Z shear b": "magenta"}
name_map = {"tensile a": "$g_{ta}$",
            "tensile b": "$g_{tb}$",
            "Z shear a": "$g_{za}$",
            "Z shear b": "$g_{zb}$"}
failure_mode_colors = np.full(shape + (4,), 0.0)
for name, gF in gFs.items():
    failure_mode_colors[minF == gF] = matplotlib.colors.to_rgba(mcd.XKCD_COLORS["xkcd:" + colormap[name]] + "c0")

FEM_color = np.asarray([0, 0, 1, .5])

def plotTwo(ax, X, Y, Z1, Z2, col1=None, col2=None):
    XX = X
    YY = Y
    ZZ = Z1
    col = col1
    if col1 is None:
        col = np.full(Z1.shape + (4,), .5, dtype=object)
        col1 = col
    if col2 is None:
        col2 = np.full(Z2.shape + (4,), .0, dtype=object)
        col2[:,:,0:4] = FEM_color

    cdict = {'red': [[0.0, 1.0, 1.0],
                     [0.5, 0.5, 0.5],
                     [1.0, 0.0, 0.0]],
             'green': [[0.0, 0.0, 0.0],
                       [0.5, 0.5, 0.5],
                       [1.0, 1.0, 1.0]],
             'blue': [[0.0, 0.0, 0.0],
                      [0.5, 0.5, 0.5],
                      [1.0, 0.0, 0.0]],
             'alpha': [[0.0, 0.0, 0.0],
                      [0.1, 0.0, 0.0],
                      [0.3, 1.0, 1.0],
                      [1.0, 1.0, 1.0]]}
    newcmp = matplotlib.colors.LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)

    cmap = plt.get_cmap('Greys')
    zz_scaled = (ZZ - np.min(ZZ)) / (np.max(ZZ) - np.min(ZZ))
    mapped = cmap(.7 - .6 * zz_scaled * zz_scaled) if compare_to_FEM else cmap(.9 - .8 * np.power(zz_scaled, 3))
    is_grey = mapped[:,:,0]
    is_grey = is_grey * 2 - 1
    is_grey *= is_grey
    mult = np.repeat(is_grey, 4, axis=1).reshape(col.shape)
    col = mapped * mult + (1 - mult) * col
    col[ZZ > .1, 3] = .6 if compare_to_FEM else 1.
    col[ZZ <= .1, 3] = 0

    if compare_to_FEM:
        col[-1, :] = 0
        XX = np.append(X, X, axis=0)
        YY = np.append(Y, Y, axis=0)
        ZZ = np.append(Z1, Z2, axis=0)
        col = np.append(col, col2, axis=0)

    ax.plot_surface(XX, YY, ZZ, facecolors=col, edgecolor='none', linewidth=0, shade=compare_to_FEM)
    #ax.scatter(X, Y, Z1, c=col.reshape(-1, 4), s=1)
    ax.set_zlim(1, 8)
    if len(wbs) < 20:
        ax.set_xticks(wbs)
    if len(Ls) < 20:
        ax.set_yticks(Ls)


if not only_one:
    if Ls.size == 1:
        fig = plt.figure()
        ax = plt.axes()
        ax.plot(wb, stress)
        ax.set_xlabel("$w_b$")
        ax.set_ylabel("stress (MPa)")
    else:
        col = failure_mode_colors
        ZZ = stress
        cmap = plt.get_cmap('Greys')
        zz_scaled = (ZZ - np.min(ZZ)) / (np.max(ZZ) - np.min(ZZ))
        mapped = cmap(.9 - .8 * zz_scaled * zz_scaled) if compare_to_FEM else cmap(.9 - .8 * np.power(zz_scaled, 3))
        is_grey = mapped[:, :, 0]
        is_grey = is_grey * 2 - 1
        is_grey *= is_grey
        mult = np.repeat(is_grey, 4, axis=1).reshape(col.shape)
        col = mapped * mult + (1 - mult) * col
        col[ZZ > .01, 3] = .6 if compare_to_FEM else 1.
        col[ZZ <= .01, 3] = 0

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.set(xlabel='$w_b$', ylabel='$L$', zlabel="stress (MPa)")
        if make_comparison_plot:
            plotTwo(ax, wb, L, stress, FEM_stress, failure_mode_colors)
            legend_elements = []
            for name, constr in gFs.items():
                color = colormap[name]
                legend_elements.append(Patch(facecolor=color, edgecolor='none', label=name_map[name]))
            legend_elements.append(Patch(facecolor=FEM_color, edgecolor='none', label='FEM'))
            ax.legend(handles=legend_elements)
        else:
            if not big or compare_to_FEM:
                Ls_ = np.linspace(np.min(Ls), np.max(Ls), 200)
                wbs_ = np.linspace(np.min(wbs), np.max(wbs), 200)
                L_, wb_ = np.meshgrid(Ls_, wbs_)
                tck = interpolate.interp2d(L, wb, stress, kind='linear')
                stress_ = tck(Ls_, wbs_)
                plot = ax.plot_surface(wb_, L_, stress_, cmap=plt.get_cmap('jet'), edgecolor='none', linewidth=0)
                fig.colorbar(plot)
            else:
                ax.scatter(wb, L, stress, c=col.reshape(-1, 4), s=1)
            if len(wbs) < 20:
                ax.set_xticks(wbs)
            if len(Ls) < 20:
                ax.set_yticks(Ls)
            if not replace_by_FEM:
                legend_elements = []
                for name, constr in gFs.items():
                    color = colormap[name]
                    legend_elements.append(Patch(facecolor=color, edgecolor='none', label=name_map[name]))
                ax.legend(handles=legend_elements)
            ax.set_zlim(3, 7)
    plt.show()