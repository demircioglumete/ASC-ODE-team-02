import numpy as np
import matplotlib.pyplot as plt

# --- files for each method ---------------------------------------------
files_cn = {
    "steps = 10":  "cn_t10.txt",
    "steps = 50":  "cn_t50.txt",
    "steps = 100": "cn_t100.txt",
}

files_imp = {
    "steps = 10":  "imp_t10.txt",
    "steps = 50":  "imp_t50.txt",
    "steps = 100": "imp_t100.txt",
}

files_impr = {
    "steps = 10":  "impr_t10.txt",
    "steps = 50":  "impr_t50.txt",
    "steps = 100": "impr_t100.txt",
}


# ---------------------------------------------------------------------
# 1) Time evolution: position vs time (Crank–Nicolson, different steps)
# ---------------------------------------------------------------------
plt.figure()
for label, fname in files_cn.items():
    data = np.loadtxt(fname, usecols=(0, 1, 2))
    t = data[:, 0]
    x = data[:, 1]
    plt.plot(t, x, label=label)

plt.xlabel("time")
plt.ylabel("position")
plt.title("Mass-Spring (Crank–Nicolson): Position for Different Time Steps")
plt.legend()
plt.grid(True)



# ---------------------------------------------------------------------
# 2) Time evolution: position vs time (Implicit, different steps)
# ---------------------------------------------------------------------
plt.figure()
for label, fname in files_imp.items():
    data = np.loadtxt(fname, usecols=(0, 1, 2))
    t = data[:, 0]
    x = data[:, 1]
    plt.plot(t, x, label=label)

plt.xlabel("time")
plt.ylabel("position")
plt.title("Mass-Spring (Implicit Euler Method): Position for Different Time Steps")
plt.legend()
plt.grid(True)


# ---------------------------------------------------------------------
# 3) Time evolution: position vs time (Improved, different steps)
# ---------------------------------------------------------------------
plt.figure()
for label, fname in files_impr.items():
    data = np.loadtxt(fname, usecols=(0, 1, 2))
    t = data[:, 0]
    x = data[:, 1]
    plt.plot(t, x, label=label)

plt.xlabel("time")
plt.ylabel("position")
plt.title("Mass-Spring (Improved Euler Method): Position for Different Time Steps")
plt.legend()
plt.grid(True)


# ---------------------------------------------------------------------
# 4) Phase plot: Crank–Nicolson
# ---------------------------------------------------------------------
plt.figure()
for label, fname in files_cn.items():
    data = np.loadtxt(fname, usecols=(0, 1, 2))
    x = data[:, 1]
    v = data[:, 2]
    plt.plot(x, v, label=label)

plt.xlabel("position")
plt.ylabel("velocity")
plt.title("Mass-Spring (Crank–Nicolson): Phase Plot for Different Time Steps")
plt.legend()
plt.grid(True)


# ---------------------------------------------------------------------
# 5) Phase plot: Implicit Euler
# ---------------------------------------------------------------------
plt.figure()
for label, fname in files_imp.items():
    data = np.loadtxt(fname, usecols=(0, 1, 2))
    x = data[:, 1]
    v = data[:, 2]
    plt.plot(x, v, label=label)

plt.xlabel("position")
plt.ylabel("velocity")
plt.title("Mass-Spring (Implicit Euler): Phase Plot for Different Time Steps")
plt.legend()
plt.grid(True)


# ---------------------------------------------------------------------
# 6) Phase plot: Improved Euler
# ---------------------------------------------------------------------
plt.figure()
for label, fname in files_impr.items():
    data = np.loadtxt(fname, usecols=(0, 1, 2))
    x = data[:, 1]
    v = data[:, 2]
    plt.plot(x, v, label=label)

plt.xlabel("position")
plt.ylabel("velocity")
plt.title("Mass-Spring (Improved Euler): Phase Plot for Different Time Steps")
plt.legend()
plt.grid(True)

plt.show()
