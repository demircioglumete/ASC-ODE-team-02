import numpy as np
import matplotlib.pyplot as plt

# Load your RK4 output
data = np.loadtxt("rk4_output.txt")

t = data[:, 0]
y_num = data[:, 1]
z_num = data[:, 2]

# Exact solution
y_exact = np.cos(t)
z_exact = -np.sin(t)

# Plot y(t)
plt.figure(figsize=(10,5))
plt.plot(t, y_num, label="RK4 y(t)", linewidth=2)
plt.plot(t, y_exact, "--", label="Exact cos(t)")
plt.xlabel("t")
plt.ylabel("y")
plt.legend()
plt.title("RK4 vs exact solution for y(t)")
plt.grid()

plt.show()
