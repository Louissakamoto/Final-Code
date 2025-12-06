'''
This file creates internal shear, bending moment and torque diagrams for a wing half-span. The wing is modelled as a rigidly mounted beam subject to 3 dimensional loading
'''
'''
Assumptions:
-wing weight is a function of the chord squared 
-engine and fuselage weights are point forces
-fuselage connected at start of halfwing span
-no torsion due to wing structures

It has been assumed that the weight distribution is W(y) = Kc^2(y), makes sense as K is in [N/m^2] and the area A = c(y)t(y), 
where t(y) = 0.14c(y) (t/c = 0.14). So we need to find the value of the constant K by integrating W(y) and saying it equals to
the total wing weight!
'''

import scipy as sp
import numpy as np
import AerodynamicLoads as al

# Parameters:
g=9.81  #[m/s^2]
WL=6030  #[N/m^2]
b=44.3  #[m]
engine_mass=118197*(0.101 + 0.023)/2  #[kg]
engine_thrust = 233591  #[N]
start_wet_wing=0
end_wet_wing=0.7
fuel_mass=41161/2  #[kg]
wing_mass=0.142*118197/2 #[kg]
front_spar=0.2*b/2 #m  as it is 20 percent of the half wingspan
rear_spar=0.65*b/2 #m  as it is 65 percent of the half wingspan
fuel_density=807.5 #kg/m^3
c_r=6.53 #[m]
c_t=1.91 #[m]
ct=0.14 #thickness to chord ratio

#Wing structural weight distribution
def il_wing_struc_dist(wsd_y, wsd_span, wsd_wing_mass, wsd_g):
    wsd_start_load = 2*wsd_wing_mass/(wsd_span/2)*wsd_g  #Assuming a linear distribution ending at zero, use the area of a triangle to find the initial load
    return -wsd_start_load/(wsd_span/2)*wsd_y + wsd_start_load  #Distribution as a function of y

#Fuel weight distribution
def il_wing_fuel_dist(fuel_y, t_to_c, root, tip, fuel_b, rho, fuel_g ):
    k=0.75  #Initial wing box factor, area of wing box to area of airfoil, can calculate more in debth but for now this should suffice
    if fuel_y <= fuel_b:
        return k*t_to_c*(root-(root-tip)*fuel_y/(fuel_b/2))**2*rho*fuel_g  #Distribution as function of y
    else:
        return 0
fuel_estimate, fuel_error = sp.integrate.quad(il_wing_fuel_dist, front_spar, rear_spar, args=(ct, c_r, c_t, b, fuel_density,g ))