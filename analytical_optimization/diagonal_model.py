import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib._color_data as mcd
import matplotlib.colors

from math import sqrt

big = True

Ls = np.linspace(.6, 3.6, 100)
wbs = np.linspace(.6, 4.8, 100)
if not big:
    Ls = 3.6
    wbs = np.linspace(.0, 4.8, 6)

L, wb = np.meshgrid(Ls, wbs)


wmin = 0.3
r = .15
wa = 2 * wmin
h = .2

sa = 47
sb = 10.5
saz = 33
sbz = 10.6

w = wa + wb


M = L - 2*r
div = np.sqrt(np.maximum(0.0001, 4*M**2 - w**2))
div[np.logical_not(np.isfinite(div))] = 999999
d = 2 * M * w / div

div = np.sqrt((w/d)**2 + 3*(d/2/M)**2)
div[np.logical_not(np.isfinite(div))] = 99999
F_tpla = sa * wa * h / div
F_pp = sb * wb * d / w
F_pp2 = sb * wb * h / div

if Ls.size == 1:
    print(F_tpla)
    print(F_pp2)

F = np.minimum.reduce([F_tpla, F_pp, F_pp2])

np.where(np.isnan(F), F, 0)

stress = F / d / w

best_wb_idx = np.unravel_index(np.argmax(stress), stress.shape)
print(f"max stress: {stress.max()} at wb: {wb[best_wb_idx]}")

if Ls.size == 1:
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(wb, stress)
    ax.set_xlabel("$w_b$")
    ax.set_ylabel("stress (MPa)")
else:
    col = np.zeros(wb.shape + (4,))
    col[:, :, 3] = 1
    col[F == F_tpla, 0] = 1
    col[F == F_pp, 1] = 1
    col[F == F_pp2, 2] = 1
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set(xlabel='$w_b$', ylabel='$L$', zlabel="stress (MPa)")
    #ax.plot_surface(wb, L, stress)
    ax.scatter(wb, L, stress, c=col.reshape(-1, 4), s=1)
plt.show()