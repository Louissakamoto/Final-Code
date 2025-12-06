from math import cos, sin, sqrt
import matplotlib.pyplot as plt
import var
import fn

y = 0

x_string = []
y_string = []

centroid, p1, p2, p3, p4 = fn.WP4_2_wingbox_shape(y)

x_corners = [p1[0], p2[0], p3[0], p4[0]]
y_corners = [p1[1], p2[1], p3[1], p4[1]]

y_step_lower = abs((p2[1]-p1[1])/(var.n_string_lower-1))
x_step_lower = abs((p2[0]-p1[0])/(var.n_string_lower-1))
for n in range(var.n_string_lower):
    x_string.append(p1[0] + n*x_step_lower)
    y_string.append(p1[1] + n*y_step_lower)

y_step_upper = abs((p4[1]-p3[1])/(var.n_string_upper-1))
x_step_upper = abs((p3[0]-p4[0])/(var.n_string_upper-1))
for n in range(var.n_string_upper):
    x_string.append(p4[0] + n*x_step_upper)
    y_string.append(p4[1] - n*y_step_upper)

print(centroid)

plt.plot(x_corners, y_corners, color='darkblue')
plt.plot(x_corners + [x_corners[0]], y_corners + [y_corners[0]], color='darkblue')
plt.plot(x_string, y_string, 'o', color='orange', label='stringers')
plt.plot(centroid[0], centroid[1], 'o', color='darkred', label='centroid')
plt.title(f'Wignbox Geometry at {y} [m]')
plt.legend()
plt.xlabel('[m]')
plt.ylabel('[m]')
plt.show()