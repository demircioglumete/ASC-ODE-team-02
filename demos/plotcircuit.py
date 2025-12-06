import numpy as np
import matplotlib.pyplot as plt

# Name of the data file produced by test_ode_circuit

datafile = "output_test_ode_circuit.txt"

# Load columns: t, U_C, U_0
data = np.loadtxt(datafile, usecols=(0, 1, 2))

t  = data[:, 0]
Uc = data[:, 1]
U0 = data[:, 2]

# ----------------------------------------------------------------------
# 1) Plot U_C(t) and U_0(t) over time
# ----------------------------------------------------------------------
plt.figure()
plt.plot(t, Uc, label='U_C(t)  (capacitor)')
plt.plot(t, U0, label='U_0(t)  (source)', linestyle='--')
plt.xlabel('time t')
plt.ylabel('voltage')
plt.title('RC Circuit: Capacitor and Source Voltage vs Time')
plt.legend()
plt.grid(True)

# ----------------------------------------------------------------------
# 2) Optional: Phase-like plot U_C vs U_0
# ----------------------------------------------------------------------
plt.figure()
plt.plot(U0, Uc)
plt.xlabel('U_0(t)  (source)')
plt.ylabel('U_C(t)  (capacitor)')
plt.title('RC Circuit: U_C vs U_0')
plt.grid(True)

plt.show()
