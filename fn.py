from math import cos, sin, sqrt, acos
import var
import numpy as np
import scipy as sp

def WP4_2_wingbox_shape(y):
    theta = 0.0497367 # in radians
    c = 6.53 - (4.62/var.b/2)*y # chord length
    d = 0.45*c # distance between spars
    s_f = 0.128*c # length of the front spar
    s_b = 0.1*c # length of the back spar

    # positions of the corner points
    p1 = (0, 0)
    p2 = (cos(theta)*d, sin(theta)*d)
    p3 = (cos(theta)*d, sin(theta)*d + s_b)
    p4 = (0, s_f)

    # y positions of the centroids for each segment
    y1 = (p1[1]+p2[1])/2
    y2 = (p2[1]+p3[1])/2
    y3 = (p3[1]+p4[1])/2
    y4 = p4[1]/2

    # x positions of the centroids for each segment
    x1 = d/2
    x2 = d
    x3 = d/2
    x4 = 0

    # lengths of each segment
    l1 = p2[1]/sin(theta)
    l2 = s_b
    l3 = sqrt((p4[0]-p3[0])**2 + (p4[1]-p3[1])**2)
    l4 = s_f

    # areas of each segment
    A1 = var.t_skin*l1
    A2 = var.t_spar*l2
    A3 = var.t_skin*l3
    A4 = var.t_spar*l4

    A_sides = A1 + A2 + A3 + A4

    # contribution of stringers to centroid
    y_step_lower = abs((p2[1]-p1[1])/(var.n_string_lower-1))
    x_step_lower = abs((p2[0]-p1[0])/(var.n_string_lower-1))
    A_y_lower = 0
    A_x_lower = 0
    for n in range(var.n_string_lower):
        A_y_lower += (p1[1] + n*y_step_lower)*var.A_stringer
        A_x_lower += (p1[0] + n*x_step_lower)*var.A_stringer

    y_step_upper = abs((p4[1]-p3[1])/(var.n_string_upper-1))
    x_step_upper = abs((p3[0]-p4[0])/(var.n_string_upper-1))
    A_y_upper = 0
    A_x_upper = 0
    for n in range(var.n_string_upper):
        A_y_upper += (p4[1] - n*y_step_upper)*var.A_stringer
        A_x_upper += (p4[0] + n*x_step_upper)*var.A_stringer
        
    A_stringers = var.A_stringer*(var.n_string_lower+var.n_string_upper)
    y_bar = (A1*y1 + A2*y2 + A3*y3 + A4*y4 + A_y_lower + A_y_upper)/(A_sides + A_stringers)
    x_bar = (A1*x1 + A2*x2 + A3*x3 + A4*x4 + A_x_lower + A_x_upper)/(A_sides + A_stringers)

    centroid = (x_bar, y_bar)

    return [centroid, p1, p2, p3, p4]
    
def WP4_2_Torsional_Stiffness(y):
    list1=WP4_2_wingbox_shape(y) #finds the properties of the cross section
    fs=list1[4][1]-list1[1][1] #height of front spar
    rs=list1[3][1]-list1[2][1] #height of rear spar
    L_top=sqrt((list1[4][0]-list1[3][0])**2+(list1[4][1]-list1[3][1])**2)
    L_bottom=sqrt((list1[2][0]-list1[1][0])**2+(list1[2][1]-list1[1][1])**2)

    h_trapezoid=list1[2][0]-list1[1][0] 

    A=(fs+rs)*h_trapezoid/2 #Cross sectional area

    integral=(fs+rs)/var.t_spar+(L_top+L_bottom)/var.t_skin #Calculates the integral in the formula for J

    J=4*A**2/integral #Polar moment of inertia

    torsional_stiffness=var.G*J
    return torsional_stiffness

def WP4_2_Ixx(y):
    list1=WP4_2_wingbox_shape(y)  
    Ixx=0

    #Calculating cross-sectional properties
    fs=list1[4][1]-list1[1][1] #height of front spar
    rs=list1[3][1]-list1[2][1] #height of rear spar
    L_top=sqrt((list1[4][0]-list1[3][0])**2+(list1[4][1]-list1[3][1])**2)
    L_bottom=sqrt((list1[2][0]-list1[1][0])**2+(list1[2][1]-list1[1][1])**2)
    h_trapezoid=list1[2][0]-list1[1][0]
    theta_top=acos(h_trapezoid/L_top)
    theta_bottom=acos(h_trapezoid/L_bottom)

    #Calculating Ixx for the front spar and rear spar
    Ixx+=fs**3*var.t_spar/12+fs*var.t_spar*((list1[4][1]+list1[1][1])/2-list1[0][1])**2 #Effect on the front spar
    Ixx+=rs**3*var.t_spar/12+rs*var.t_spar*((list1[3][1]+list1[2][1])/2-list1[0][1])**2 #Effect on the rear spar

    #Calculating Ixx for the top and bottom segments
    Ixx+=var.t_skin*L_top**3*sin(theta_top)**2/12+var.t_skin*L_top*((list1[4][1]+list1[3][1])/2-list1[0][1])**2 #Top segment
    Ixx+=var.t_skin*L_bottom**3*sin(theta_bottom)**2/12+var.t_skin*L_bottom*((list1[2][1]+list1[1][1])/2-list1[0][1])**2 #Bottom segment

    #Calculate the effect of the stringers
    lower_difference=abs(list1[2][1]-list1[1][1])/(var.n_string_lower-1)
    for i in range(var.n_string_lower):
       z_lower=list1[1][1]+lower_difference*i #Add because it is going up (left to right)
       Ixx+=var.A_stringer*(list1[0][1]-z_lower)**2
    
    upper_difference=abs(list1[4][1]-list1[3][1])/(var.n_string_upper-1)
    for i in range(var.n_string_upper):
        z_upper=list1[4][1]-upper_difference*i #Subtract because it is going down (left to right)
        Ixx+=var.A_stringer*(list1[0][1]-z_upper)**2

    return Ixx

def WP4_2_Izz(y):
    list1=WP4_2_wingbox_shape(y)  
    Izz=0

    #Calculating cross-sectional properties
    fs=list1[4][1]-list1[1][1] #height of front spar
    rs=list1[3][1]-list1[2][1] #height of rear spar
    L_top=sqrt((list1[4][0]-list1[3][0])**2+(list1[4][1]-list1[3][1])**2)
    L_bottom=sqrt((list1[2][0]-list1[1][0])**2+(list1[2][1]-list1[1][1])**2)
    h_trapezoid=list1[2][0]-list1[1][0]
    theta_top=acos(h_trapezoid/L_top)
    theta_bottom=acos(h_trapezoid/L_bottom)

    #Calculating Ixx for the front spar and rear spar
    Izz+=fs*var.t_spar*(list1[4][0]-list1[0][0])**2 #Effect on the front spar
    Izz+=rs*var.t_spar*(list1[3][0]-list1[0][0])**2 #Effect on the rear spar

    #Calculating Ixx for the top and bottom segments
    Izz+=var.t_skin*L_top**3*cos(theta_top)**2/12+var.t_skin*L_top*((list1[4][0]+list1[3][0])/2-list1[0][0])**2 #Top segment
    Izz+=var.t_skin*L_bottom**3*cos(theta_bottom)**2/12+var.t_skin*L_bottom*((list1[2][0]+list1[1][0])/2-list1[0][0])**2 #Bottom segment

    #Calculate the effect of the stringers
    #Calculate the effect of the stringers
    lower_difference=abs(list1[2][0]-list1[1][0])/(var.n_string_lower-1)
    for i in range(var.n_string_lower):
       x_lower=list1[1][0]+lower_difference*i #Add because it is going right
       Izz+=var.A_stringer*(list1[0][0]-x_lower)**2
    
    upper_difference=abs(list1[4][0]-list1[3][0])/(var.n_string_upper-1)
    for i in range(var.n_string_upper):
        x_upper=list1[4][0]+upper_difference*i #Add because it is going right
        Izz+=var.A_stringer*(list1[0][0]-x_upper)**2

    return Izz
