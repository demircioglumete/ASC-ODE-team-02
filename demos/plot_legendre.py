import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("legendre_data.txt", skiprows=1)

x = data[:,0]

# P0..P5 and derivatives
P = [data[:, 1 + 2*n] for n in range(6)]
dP = [data[:, 2 + 2*n] for n in range(6)]

# ------- Plot Polynomials -------
plt.figure(figsize=(10,6))
for n in range(6):
    plt.plot(x, P[n], label=f"P{n}(x)")
plt.title("Legendre Polynomials P0–P5")
plt.xlabel("x")
plt.ylabel("P_n(x)")
plt.grid(True)
plt.legend()
plt.show()

# ------- Plot Derivatives -------
plt.figure(figsize=(10,6))
for n in range(6):
    plt.plot(x, dP[n], label=f"P{n}'(x)")
plt.title("Derivatives of Legendre Polynomials")
plt.xlabel("x")
plt.ylabel("P_n'(x)")
plt.grid(True)
plt.legend()
plt.show()
