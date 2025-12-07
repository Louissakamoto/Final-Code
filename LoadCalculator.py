import numpy as np
import scipy as sp

#Finding the structual weight distribution as a function of y
def weightdistribution2(mass, span, ch, g):
    def structural_integrand(y):
        return ch(y)**2 
    
    integral_result, _ = sp.integrate.quad(structural_integrand, 0, span / 2)
    
    if isinstance(integral_result, (tuple, list)):
        integral_value = integral_result[0]
    
    else:
        integral_value = integral_result
    K = mass / integral_value
    def Wd(y):
        return K * structural_integrand(y) * g
    
    return Wd

#Fuel weight distribution
def il_wing_fuel_dist(fuel_y, wet_span, root, tip, fuel_b, rho, fuel_g, bol ):
    if bol:
        if fuel_y <= wet_span:
            return rho * fuel_g * 0.051525*(root + (tip - root)*fuel_y*2/fuel_b)**2  #Distribution as function of y
        else:
            return 0
    else:
        return 0

#Finding engine weight
def engine_weight(ew_engine_mass, ew_g, ew_y_engine):
    return (ew_y_engine, ew_engine_mass*ew_g)  #Returns tuple with spanwise position of force and magnitude

#A functoin that computes the resultant load
def load(span, g, c_to_t, root_chord, tip_chord, end_wet, fuel_rho, angle_attack, wed, led, ded, check, sf):

    def structualfixed(y):
        return wed(y)*np.cos(angle_attack)
    
    def fuelfixed(y):
        return il_wing_fuel_dist(y, end_wet, root_chord, tip_chord, span, fuel_rho, g, check)*np.cos(angle_attack)

    def combined(y):
        return sf* (-(led(y)*np.cos(angle_attack)+ded(y)*np.sin(angle_attack)) + structualfixed(y) + fuelfixed(y))
    
    return combined


