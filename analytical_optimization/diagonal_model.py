import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib._color_data as mcd
import matplotlib.colors

from math import sqrt
from math import pi

big = False
only_one = False

compare_to_FEM = False

Ls = np.linspace(.6, 3.6, 200)
wbs = np.linspace(.6, 4.8, 200)
if not big:
    Ls = np.asarray([3.6])
    wbs = np.linspace(.6, 3.6, 6)
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
        #ax.plot_surface(wb, L, stress)
        ax.scatter(wb, L, stress, c=col.reshape(-1, 4), s=1)
        legend_elements = []
        for name, constr in gFs.items():
            color = colormap[name]
            legend_elements.append(Patch(facecolor=color, edgecolor='none', label=name_map[name]))
        ax.legend(handles=legend_elements)
        ax.set_zlim(0, 7)
    plt.show()