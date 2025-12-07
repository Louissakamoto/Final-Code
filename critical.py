import Load_distribution as LD
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

load_cases = [

    # For the V-n diagrams, refer to appendix A,
    #The first true statement means fuel is included, the second true statement means max payload included.

    # a
    #V_eas, altitude, n, fuel?, payload?
    #in reality the altitude is not required because the dynamic pressure is computed using V_EAS, but its included anyways.
    (163, 9448, 2.5, True, True),   #0
    (154, 9448, -1.0, True, True),    #1
    (152, 9448, 2.5, True, True),   #2

    # b
    (193, 0, 2.5, False, True),  #3
    (154, 0, -1.0, False, True),   #4
    (123, 0, 2.5, False, True),   #5

    # c
    (193, 0, 2.5, False, False),   #6
    (154, 0, -1.0, False, False),   #7
    (107, 0, 2.5, False, False),   #8
    
    #d
    (163, 9448, 2.5, False, True),   #9
    (154, 9448, -1, False, True),   #10
    (123, 9448, 2.5, False, True),    #11
    
    #e
    (163, 9448, 2.5, False, False),   #12
    (154, 9448, -1, False, False),    #13
    (107, 9448, 2.5, False, False),    #14

    #f
    (193, 0, 2.5, True, True),    #15
    (154, 0, -1.0, True, True),    #16
    (152, 0, 2.5, True, True),    #17
    
]

maxS = [0, 0]
maxM = [0, 0]
maxT = [0, 0]
countS = [0, 0]
countM = [0, 0]
countT = [0, 0]
count = 0

plt.figure(figsize=(12, 8))
plt.subplots_adjust(wspace=0.5)
plt.subplots_adjust(hspace=0.5)
plt.suptitle("Critical Loading Analysis")

for case in load_cases:
    V_eas, alt, n, fuel, payload = case
    x, S, M, T= LD.internalcalculator(V_eas, alt, n, fuel, payload)

    if n > 0:
        if np.max(np.abs(S)).item() > maxS[0]:
            maxS[0] = np.max(np.abs(S))
            countS[0] = count
        if np.max(np.abs(M)).item() > maxM[0]:
            maxM[0] = np.max(np.abs(M))
            countM[0] = count
        if np.max(np.abs(T)).item() > maxT[0]:
            maxT[0] = np.max(np.abs(T))
            countT[0] = count

        plt.subplot(231)
        plt.plot(x, S)
        plt.title("Shear, positive load factor")
        plt.xlabel("Spanwise position y [m]")
        plt.ylabel("Shear force S(y) [N]")
        

        plt.subplot(232)
        plt.plot(x, M)
        plt.title("Bending Moment, positive load factor")
        plt.xlabel("Spanwise position y [m]")
        plt.ylabel("Moment M(y) [Nm]")

        plt.subplot(233)
        plt.plot(x, T)
        plt.title("Torsion, positive load factor")
        plt.xlabel("Spanwise position y [m]")
        plt.ylabel("Torsion T(y) [Nm]")

    elif n < 0:
        if np.max(np.abs(S)).item() > maxS[1]:
            maxS[1] = np.max(np.abs(S))
            countS[1] = count
        if np.max(np.abs(M)).item() > maxM[1]:
            maxM[1] = np.max(np.abs(M))
            countM[1] = count
        if np.max(np.abs(T)).item() > maxT[1]:
            maxT[1] = np.max(np.abs(T))
            countT[1] = count

        plt.subplot(234)
        plt.plot(x, S)
        plt.title("Shear, negative load factor")
        plt.xlabel("Spanwise position y [m]")
        plt.ylabel("Shear force S(y) [N]")

        plt.subplot(235)
        plt.plot(x, M)
        plt.title("Bending Moment, negative load factor")
        plt.xlabel("Spanwise position y [m]")
        plt.ylabel("Moment M(y) [Nm]")

        plt.subplot(236)
        plt.plot(x, T)
        plt.title("Torsion, negative load factor")
        plt.xlabel("Spanwise position y [m]")
        plt.ylabel("Torsion T(y) [Nm]")

    count += 1
    print("Im looping :) ", count,"/, please wait")


plt.subplot(231)
plt.gca().invert_yaxis()
plt.subplot(234)
plt.gca().invert_yaxis()

plt.show()

print(maxS, maxM, maxT, countS, countM, countT)

#positive loading
V_easp, altp, npos, fuelp, payloadp = load_cases[15]
x, Sp, Mdummy, Tdummy = LD.internalcalculator(V_easp, altp, npos, fuelp, payloadp)
V_easp, altp, npos, fuelp, payloadp = load_cases[15]
x, Sdummy, Mp, Tdummy = LD.internalcalculator(V_easp, altp, npos, fuelp, payloadp)
V_easp, altp, npos, fuelp, payloadp = load_cases[2]
x, Sdummy, Mdummy, Tp = LD.internalcalculator(V_easp, altp, npos, fuelp, payloadp)

Spf = sp.interpolate.interp1d(x, Sp, kind="cubic", fill_value='extrapolate')
Mpf = sp.interpolate.interp1d(x, Mp, kind="cubic", fill_value='extrapolate')
Tpf = sp.interpolate.interp1d(x, Tp, kind="cubic", fill_value='extrapolate')

#negative loading
V_easp, altp, npos, fuelp, payloadp = load_cases[1]
x, Sn, Mdummy, Tdummy = LD.internalcalculator(V_easp, altp, npos, fuelp, payloadp)
V_easp, altp, npos, fuelp, payloadp = load_cases[1]
x, Sdummy, Mn, Tdummy = LD.internalcalculator(V_easp, altp, npos, fuelp, payloadp)
V_easp, altp, npos, fuelp, payloadp = load_cases[1]
x, Sdummy, Mdummy, Tn = LD.internalcalculator(V_easp, altp, npos, fuelp, payloadp)

Snf = sp.interpolate.interp1d(x, Sn, kind="cubic", fill_value='extrapolate')
Mnf = sp.interpolate.interp1d(x, Mn, kind="cubic", fill_value='extrapolate')
Tnf = sp.interpolate.interp1d(x, Tn, kind="cubic", fill_value='extrapolate')


#plot critical shear
plt.figure(figsize=(13, 5))
plt.subplots_adjust(wspace=0.5)
plt.suptitle("Shear, Critical Load Cases")

plt.subplot(121)
plt.plot(x, Sp)
plt.title("Positive Load Case")
plt.xlabel("Span Wise Location [m]")
plt.ylabel("Shear [N]")
plt.gca().invert_yaxis()

plt.subplot(122)
plt.plot(x, Sn)
plt.title("Negative Load Case")
plt.xlabel("Span Wise Location [m]")
plt.ylabel("Shear [N]")
plt.gca().invert_yaxis()

plt.show()

#plot critical bending
plt.figure(figsize=(13, 5))
plt.subplots_adjust(wspace=0.5)
plt.suptitle("Bending Moment, Critical Load Cases")

plt.subplot(121)
plt.plot(x, Mp)
plt.title("Positive Load Case")
plt.xlabel("Span Wise Location [m]")
plt.ylabel("Bending Moment [Nm]")

plt.subplot(122)
plt.plot(x, Mn)
plt.title("Negative Load Case")
plt.xlabel("Span Wise Location [m]")
plt.ylabel("Bending Moment [Nm]")

plt.show()

#plot critical torsion
plt.figure(figsize=(13, 5))
plt.subplots_adjust(wspace=0.5)
plt.suptitle("Torque, Critical Load Cases")

plt.subplot(121)
plt.plot(x, Tp)
plt.title("Positive Load Case")
plt.xlabel("Span Wise Location [m]")
plt.ylabel("Torque [Nm]")

plt.subplot(122)
plt.plot(x, Tn)
plt.title("Negative Load Case")
plt.xlabel("Span Wise Location [m]")
plt.ylabel("Torque [Nm]")

plt.show()

