import math

def ISA(h):
    g = 9.80665
    R = 287
    T0 = 288.15
    p0 = 101325.0
    a = -0.0065
    rho0 = p0/(R * T0)

    layers = [
        (0, 11000, -0.0065),
        (11000, 20000, 0.0),
        (20000, 32000, 0.0010),
        (32000, 47000, 0.0028),
        (47000, 51000, 0.0),
        (51000, 71000, -0.0028),
        (71000, 86000, -0.0020)
    ]

    T = T0
    p = p0
    rho = rho0

    for h_base, h_top, lapse in layers:
        if h <= h_top:
            delta = h - h_base
            if lapse == 0:
                p = p * math.exp(-g * delta / (R * T))
            else:
                T_new = T + lapse * delta
                p = p * (T_new / T) ** (-g / (lapse * R))
                T = T_new
            rho = p / (R * T)
            break
        else:
            delta = h_top - h_base
            if lapse == 0:
                p = p * math.exp(-g * delta / (R * T))
            else:
                T_new = T + lapse * delta
                p = p * (T_new / T) ** (-g / (lapse * R))
                T = T_new

    return rho
