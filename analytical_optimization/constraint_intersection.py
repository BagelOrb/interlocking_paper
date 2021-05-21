
from math import *

def debug(variable):
    print(variable, '=', repr(eval(variable)))

sa = 47
sb = 10.5
saz = 33
sbz = 9.0

ta = .5 * sa
tb = .5 * sb
taz = .5 * saz
tbz = .5 * sbz


line_w = .3
layer_thickness = .1;

wa_min = 2 * line_w
wb_min = 2 * line_w
w_max = 6 * 2 * wa_min
l_max = 6 * 2 * wa_min

h_min = 2 * layer_thickness
h_max = 6 * h_min

wa = wb = va = vb = hc = hf = F = inf

def params1():
    global wa, wb, va, vb, hc, hf, F
    det = 1 + 4*l_max / wa_min

    F1 = sb*h_min / 3 * (  2*l_max + wa_min * ( 1 + sqrt(det) )  )
    F2 = sb*h_min / 3 * (  2*l_max + wa_min * ( 1 - sqrt(det) )  )

    # print(F1)
    # print(F2)

    F = F2

    vb = 3*F / (4 * tb*h_min)
    wb = 2*wa_min*sa/sb
    hf = F / (sa * 2 * wa_min)
    va = l_max - vb
    hc = h_min
    wa = 2* wa_min


def params2():
    global wa, wb, va, vb, hc, hf, F

    vb = -.5 * sqrt(wa_min*(wa_min+4*l_max)) + l_max + .5 * wa_min
    hc = 3 * wa_min * h_min * sa / vb / sb
    wb = 2 * wa_min * sa / sb
    F = 2 * wa_min * h_min * sa
    va = l_max - vb
    hf = h_min
    wa = 2 * wa_min

def fem_params(rz = 3, rw = 2, rx = 3.5, ry = 1.0, FF = 5.77):
    global wa, wb, va, vb, hc, hf, F
    F = FF
    wa = 2 * wa_min
    hc = h_min
    wb = wa * rx
    va = wa * rw
    hf = hc * rz
    vb = va * ry

#fem_params( rz= 1.8106122448979591, rw= 1.0956410256410258, rx= 1.8056179775280898, ry= 4.475842696629214, FF= 4.119285264847512)
fem_params()

print(f'{wa = }; {wb = }; {va = }; {vb = }; {hf = }; {hc = }; {F = };')
print(f' rz={hf/hc = }, rw={va/wa = }, rx={wb/wa = }, ry={vb/va = }')

g1 = 1 - wa / 2 / wa_min
g2 = 1 - wb / 2 / wb_min
g3 = 1 - va / wa_min
g4 = 1 - vb / wb_min
g5 = 1 - hf / h_min
g6 = 1 - hc / h_min
g7 = (va + vb) / l_max - 1
g8 = F / (wa * hf * sa) - 1  # tensile
g9 = F / (wb * hf * sb) - 1
g10 = 3 * F / (4 * va * hc * ta) - 1  # shear
g11 = 3 * F / (4 * vb * hc * tb) - 1
g12 = 3 * F / (4 * va * wa * taz) - 1  # z shear
g13 = 3 * F / (4 * vb * wa * tbz) - 1
g14 = 3 * F * wb / (4 * va * va * hc * sa) - 1  # bending
g15 = 3 * F * wa / (4 * vb * vb * hc * sb) - 1

print(f'{g1 = }, {g2 = }, {g3 = }, {g4 = }, {g5 = }, {g6 = }\n, {g7 = }, {g8 = }, {g9 = }, {g10 = }, {g11 = }, {g12 = }, {g13 = }, {g14 = }, {g15 = }')

objective = F / ( (wa + wb) * (hf + hc) )

print("Total effective stress: ", objective)

print(sqrt(h_min) / sqrt(3*wa*(sa+sb)/(l_max*h_min*sb*sa)) / sqrt(sb) )