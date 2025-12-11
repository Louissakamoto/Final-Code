from math import cos, sin, sqrt, acos, tan, atan
import var
import numpy as np
import matplotlib.pyplot as plt

def WP4_2_wingbox_shape(y):
    theta = 0.0497367 # in radians
    c = 6.53 - (4.62/var.b/2)*y # chord length
    d = 0.45*c # distance between spars
    s_f = 0.128*c # length of the front spar
    s_b = 0.1*c # length of the back spar

    # positions of the corner points
    p1 = (0, 0)
    p2 = (d, tan(theta)*d)
    p3 = (d, tan(theta)*d + s_b)
    p4 = (0, s_f)

    #calculate the angle at the upper surface
    theta2 = atan((s_f - p3[1])/d)

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
    l1 = d/cos(theta)
    l2 = s_b
    l3 = d/cos(theta2)
    l4 = s_f

    # areas of each segment
    A1 = var.t_skin*l1
    A2 = var.t_spar*l2
    A3 = var.t_skin*l3
    A4 = var.t_spar*l4

    if var.multi:
        p5 = (var.spar_location*d, var.spar_location*d*tan(theta))
        p6 = (var.spar_location, s_f - var.spar_location*d*tan(theta2))
        y5 = (p5[1]+p6[1])/2
        x5 = var.spar_location*d
        l5 = p6[1] - p5[1]
        A5 = var.t_spar *l5

    A_sides = A1 + A2 + A3 + A4
    if var.multi:
        A_sides += A5

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
    if var.multi:
        y_bar = (A1*y1 + A2*y2 + A3*y3 + A4*y4 + A5*y5 + A_y_lower + A_y_upper)/(A_sides + A_stringers)
        x_bar = (A1*x1 + A2*x2 + A3*x3 + A4*x4 + A5*x5 + A_x_lower + A_x_upper)/(A_sides + A_stringers)

    centroid = (x_bar, y_bar)
    result = [centroid, p1, p2, p3, p4]
    if var.multi:
        result = [centroid, p1, p2, p3, p4, p5, p6]

    return result
    
def WP4_2_Torsional_Stiffness(y):

    if not var.multi:
        theta = 0.0497367 # in radians
        c = 6.53 - (4.62/var.b/2)*y # chord length
        d = 0.45*c # distance between spars
        s_f = 0.128*c # length of the front spar
        s_b = 0.1*c # length of the back spar
        t1 = var.t_spar
        t2 = var.t_skin
        l2 = s_b
        l4 = s_f
        p3 = (d, d*tan(theta)+l2)
        theta2 = atan((l4 - p3[1])/d)
        l1 = d/cos(theta)
        l3 = d/cos(theta2)

        A = (l4 + l2) * d /2
        I = l4/t1 + l1/t2 + l2/t1 + l3/t2
        GJ = var.G * 4 * A**2 / I
        return GJ
    else:
        theta = 0.0497367 # in radians
        c = 6.53 - (4.62/var.b/2)*y # chord length
        d = 0.45*c # distance between spars
        s_f = 0.128*c # length of the front spar
        s_b = 0.1*c # length of the back spar
        f = var.spar_location

        t1 = var.t_spar
        t2 = var.t_skin

        l2 = s_b
        l4 = s_f
        p3 = (d, d*tan(theta)+l2)

        theta2 = atan((l4 - p3[1])/d)

        p5 = (f*d, f*d*tan(theta))
        p6 = (f*d, l4 - f*d*tan(theta2))

        l1 = d/cos(theta)
        l3 = d/cos(theta2)

        s1 = f*l1
        s2 = (1 - f)*l1
        s3 = l2
        s4 = (1 -f)*l3
        s5 = f*l3
        s6 = l4
        s7 = p6[1] - p5[1]

        A1 = (s6 + s7) * f * d/2
        A2 = (s7 + s3) * (1 -f) * d/2

        c1 = 1/(2*A1*var.G)*(s1/t2+s7/t1+s5/t2+s6/t1)
        c2 = -1/(2*A1*var.G)*(s7/t1)
        c3 = -1/(2*A2*var.G)*(s7/t1)
        c4 = 1/(2*A1*var.G)*(s2/t2+s3/t1+s4/t2+s7/t1)

        LHS = [[0, 2*A1, 2*A2], [-1, c1, c2], [-1, c3, c4]]
        RHS = [1, 0, 0]
        sol = np.linalg.solve(LHS, RHS)
        return 1/sol[0]
    
        # p3 = (d, tan(theta)*d + s_b)
        # theta2 = atan((s_f - p3[1])/d)
        # p5 = (var.spar_location*d, var.spar_location*d*tan(theta))
        # p6 = (var.spar_location*d, s_f - var.spar_location*d*tan(theta2))
        # l5 = p6[1] - p5[1]

        # #determine the areas of the cells
        # h_1 = var.spar_location*d 
        # h_2 = (1-var.spar_location)*d

        # #determine constants for multicell torsion
        # t1 = var.t_spar
        # t2 = var.t_skin

        # G = var.G
        # s1 = l4
        # s2 = var.spar_location*l1 
        # s3 = l5
        # s4 = var.spar_location*l3
        # s5 = l1 - s2
        # s6 = l2
        # s7 = l3 - s4

        # A_1 = (s1+s3)*h_1/2
        # A_2 = (s3+s6)*h_2/2

        # k1 = 1/(2*A_1*G)*(s1/t1+s2/t2+s3/t1+s4/t2)
        # k2 = -1/(2*A_1*G)*s3/t1
        # k3 = 1/(2*A_2*G)*(s5/t2+s6/t1+s7/t2+s3/t1)
        # k4 = -1/(2*A_2*G)*s3/t1

        # lh_matrix = np.array([[0, 2*A_1,2*A_2],
        #                    [-1, k4, k3],
        #                    [-1, k1, k2]])
        # rh_matrix = np.array([1, 0, 0])
        # solution = np.linalg.solve(lh_matrix, rh_matrix)
        # torsional_stiffness = 1/solution[0]
        # return torsional_stiffness
     
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

    if var.multi:
        theta = 0.0497367 # in radians
        c = 6.53 - (4.62/var.b/2)*y # chord length
        d = 0.45*c # distance between spars
        s_f = 0.128*c # length of the front spar
        s_b = 0.1*c # length of the back spar
        p3 = (d, tan(theta)*d + s_b)
        theta2 = atan((s_f - p3[1])/d)

        p5 = (var.spar_location*d, var.spar_location*d*tan(theta))
        p6 = (var.spar_location, s_f - var.spar_location*d*tan(theta2))
        y5 = (p5[1]+p6[1])/2
        l5 = p6[1] - p5[1]
        A5 = var.t_spar *l5
        Ixx += 1/12 * var.t_spar * l5**2 + A5*(y5 - list1[0][1])**2

    return Ixx

#Create a y_list with step size 0.1
# y_list = np.zeros(224)
# for i in range(len(y_list)):
#     y_list[i] = i/10

# Ixx_list = np.zeros(224)
# for i in range(len(Ixx_list)):
#     Ixx_list[i] = WP4_2_Ixx(i/10)

# GJ_list = np.zeros(224)
# for i in range(len(GJ_list)):
#     GJ_list[i] = WP4_2_Torsional_Stiffness(i/10)

# plt.figure(figsize=(13, 5))
# plt.subplots_adjust(wspace=0.5)

# plt.subplot(121)
# plt.plot(y_list, Ixx_list)
# plt.title("Moment of Inertia")
# plt.xlabel("Span-Wise location [m]")
# plt.ylabel("Ixx [m^4]")

# plt.subplot(122)
# plt.plot(y_list, GJ_list)
# plt.title("Torsional Stiffness")
# plt.xlabel("Span-Wise location [m]")
# plt.ylabel("GJ [Nm^2]")

# plt.show()