# a = PLA
# b = PP

from math import *

# E, sigma_yield, epsilon_yield
PLA = [1820, 37, 3.1]
PP = [220, 8.7, 18]
TPU = [26, 8.6, 55]
Nylon = [579, 27.8, 20]


mat_a = PLA
mat_b = PP

Ea = mat_a[0]
Eb = mat_b[0]
sa = mat_a[1]
sb = mat_b[1]
sa_shear = sa * .75
sb_shear = sb * .75

hf = .4
hc = .2
h = hc + hf

smallest_finger_width = .6

wa_min = 0 #.35
wb_min = 0 #.38

layer_height = 0.1

def NisTwo():
    N = 2
    wa1 = 1.1
    wb1 = wa1 / (sb / sa)
    wb2 = wb1 / (1 + Eb / Ea)
    w = wa1 + wb2
    wa2 = w - wb1

    print(wa1, wb1, wb2, wa2)

    Fa1 = wa1 * sa
    Fb1 = wb1 * sb
    Fa2 = wa2 * sa
    Fb2 = Fa1 - Fa2

    assert(Fa1 == Fb1)
    assert(abs(Fa2 + Fb2  -  Fa1) < .0001)

    assert(abs(Fa2 / Fb2   -   (Ea * wa2) / (Eb * wb2))  < .00001)

    sb2 = Fb2 / wb2

    print(sa, sb, sb2)

    ea1 = Fa1 / ( Ea * h * wa1 )
    ea2 = Fa2 / ( Ea * h * wa2 )
    eb1 = Fb1 / ( Eb * h * wb1 )
    eb2 = Fb2 / ( Eb * h * wb2 )
    print( ea1, ea2, eb1, eb2)

    assert( ea2   -   eb2 < .00001)

def general_case(N = 2):
    global sa_shear
    global sb_shear

    wb = [1] # temp val for wN
    for i in range(1,N):
        wb.append(wb[-1] + pow(Eb / Ea, i))
    wb.reverse()

    wb1 = wb[0]
    wa1 = wb1 * sb / sa

    w = wb[-1] + wa1

    wa = [wa1]
    for i in range(1,N):
        wa.append(w - wb[N - 1 - i])

    scaling_factor = smallest_finger_width / wa[-1]
    for i in range(N):
        wa[i] *= scaling_factor
        wb[i] *= scaling_factor
    w *= scaling_factor
    wa1 *= scaling_factor

    #hc = 2 * wa[-1]
    #hc = floor(hc / layer_height) * layer_height

    Fa1 = wa1 * hf * sa
    Fa_beam = []
    for i in range(0, N-1):
        Fa_beam.append(sa * hf * (wa[i] - wa[i+1]))
    Fa_beam.append(wa[-1] * hf * sa)

    va = []
    vb = []
    total_length = 0
    for i in range(0, N):
        va.append(max(wa_min, Fa_beam[i] / 2 / sa_shear / hc))
        vb.append(max(wb_min, Fa_beam[N-1-i] / 2 / sb_shear / hc))
        # va.append(max(0.35, 4/3 * hf/hc * wb[-1] * pow(Eb/Ea, N-i-1)))
        # vb.append(max(0.38, sa / (.75*sb) * hf/hc * wb[-1] * pow(Eb/Ea, i)))
        total_length += va[-1] + vb[-1]



    effective_stress = Fa1 / w / h
    print("max effective stress: ",  effective_stress, "MPa")
    print("total_length: ", total_length)
    print("")
    print("element_width =", w, ";")
    print("PLA_widths =", wa, ";")
    print("PP_widths =", wb, ";")
    print("hf = ", hf, ";")
    print("hc = ", hc, ";")
    print("PLA_vs =", va, ";")
    print("PP_vs =", vb, ";")
    print("")


    Fa = [Fa1]
    for i in range(1, N):
        Fa.append(wa[i] * hf * sa)
    Fb = [Fa1]
    for i in range(1, N):
        Fb.append(Fa1 - Fa[N - i])
    sb_s = []
    for i in range(N):
        sb_s.append(Fb[i] / h / wb[i])
    print("PLA forces: ", Fa)
    print("PP forces: ", Fb)
    print("PP stresses: ", sb_s)
    print("PLA beam forces: ", Fa_beam)
    sa_beam = []
    sb_beam = []
    sa_shears = []
    sb_shears = []
    for i in range(N):
        sa_beam.append(Fa_beam[i] / (va[i] * hc))
        sb_beam.append(Fa_beam[N-1-i] / (vb[i] * hc))
        sa_shears.append(Fa_beam[i] / (2 * va[i] * wa[i]))
        sb_shears.append(Fa_beam[N-1-i] / (2 * vb[i] * wb[i]))
    print("PLA beam stresses: ", sa_beam)
    print("PP beam stresses: ", sb_beam)
    print("PLA shear stresses: ", sa_shears)
    print("PP shear stresses: ", sb_shears)





#general_case(30)
#general_case(300)

#NisTwo()


if (True):
    for i in range(1,4):
        print("\n N=",i)
        general_case(i)