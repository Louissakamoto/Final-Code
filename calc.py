from scipy.integrate import quad
import fn
import var
import critical as cr
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

b = 44.3
Mp = cr.Mpf
Tp = cr.Tpf
Mn = cr.Mnf
Tn = cr.Tnf

def PMF(y):
    return -Mp(y)/(var.E*fn.WP4_2_Ixx(y))
def NMF(y):
    return -Mn(y)/(var.E*fn.WP4_2_Ixx(y))

def PTF(y):
    return Tp(y)/fn.WP4_2_Torsional_Stiffness(y)
def NTF(y):
    return Tn(y)/fn.WP4_2_Torsional_Stiffness(y)

def bending(L_end, func, n_p):

    x = np.linspace(0, L_end, n_p)
    f = np.zeros_like(x)     

    for i in range(len(x)):
        f[i] = func(x[i]) 

    F = cumulative_trapezoid(f, x, initial=0.0)
    FF = cumulative_trapezoid(F, x, initial = 0.0)

    return x, FF

def twisting(L_end, func, n_p):
    
    x = np.linspace(0, L_end, n_p)
    f = np.zeros_like(x)     

    for i in range(len(x)):
        f[i] = func(x[i]) 

    F = cumulative_trapezoid(f, x, initial = 0.0)

    return F

x, pbending = bending(b/2, PMF, 10000)
_, nbending = bending(b/2, NMF, 10000)

ptwisting = twisting(b/2, PTF, 10000)*360/np.pi
ntwisting = twisting(b/2, NTF, 10000)*360/np.pi

print("the tip deflection for the critical positive load case is ", pbending[-1], " [m]")
print("the tip rotation for the critical positive load case is ", ptwisting[-1], " [deg]")
print("the tip deflection for the critical negative load case is ", nbending[-1], " [m]")
print("the tip rotation for the critical negative load case is ", ntwisting[-1], " [deg]")

plt.figure(figsize=(13, 5))
plt.subplots_adjust(wspace=0.5)
plt.suptitle("Positive Critical Case Deflections")

plt.subplot(121)
plt.plot(x, pbending)
plt.title("vertical deflection")
plt.xlabel("Span Wise Location [m]")
plt.ylabel("Vertical Displacement [m]")
plt.gca().invert_yaxis()

plt.subplot(122)
plt.plot(x, ptwisting)
plt.title("torsional deflection")
plt.xlabel("Span Wise Location [m]")
plt.ylabel("Angular Dispalcement [deg]")

plt.show()

plt.figure(figsize=(13, 5))
plt.subplots_adjust(wspace=0.5)
plt.suptitle("Negative Critical Case Deflections")

plt.subplot(121)
plt.plot(x, nbending)
plt.title("vertical deflection")
plt.xlabel("Span Wise Location [m]")
plt.ylabel("Vertical Displacement [m]")
plt.gca().invert_yaxis()

plt.subplot(122)
plt.plot(x, ntwisting)
plt.title("torsional deflection")
plt.xlabel("Span Wise Location [m]")
plt.ylabel("Angular Dispalcement [deg]")

plt.show()

# '''method which uses sp.integrate.quad for the qurious TA'''
# def dvdy(bending):

#     def result(y):
#         return quad(bending, 0, y)[0]

#     return result

# def v(dbending):

#     def result(y):
#         return quad(dbending, 0, y)[0]

#     return result

# ddisplacement = dvdy(PMF)
# displacement = v(ddisplacement)

# displacement_x, error_x = quad(ddisplacement, 0, var.b/2)
# displacement_theta, error_theta = quad(PTF, 0, var.b/2)
