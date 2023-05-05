import numpy as np
import matplotlib.pyplot as plt

# Constants
a = 1 # m
b = 0.8 # m
mu_s = 0.262 # kg/m^2
T = 3990 # N/m

# Natural frequencies
n_max = 10
omega_nm = np.zeros((n_max, n_max))
for m in range(1, n_max+1):
    for n in range(1, n_max+1):
        omega_nm[m-1, n-1] = np.pi*np.sqrt((m/a)**2 + (n/b)**2)*np.sqrt(T/mu_s)

# Excitation function
def f(x, y):
    return np.ones_like(x)

# Spatial grid
N = 100
x = np.linspace(0, a, N)
y = np.linspace(0, b, N)
X, Y = np.meshgrid(x, y)

# Calculate modal responses
psi_nm = np.zeros((n_max, n_max))
for m in range(1, n_max+1):
    for n in range(1, n_max+1):
        psi_nm[m-1, n-1] = np.sum(np.sin(m*np.pi*X/a)*np.sin(n*np.pi*Y/b)*f(X, Y))

# Calculate transfer function
G = np.zeros((n_max, n_max), dtype=np.complex)
for m in range(1, n_max+1):
    for n in range(1, n_max+1):
        G[m-1, n-1] = psi_nm[m-1, n-1]/(mu_s*(omega_nm[m-1, n-1])**2)

# Plot results
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
fig.suptitle('Frequency response of rectangular membrane')
axs[0, 0].set_title('Excitation')
axs[0, 0].pcolormesh(X, Y, f(X, Y), shading='auto')
axs[0, 1].set_title('Modal responses')
axs[0, 1].pcolormesh(np.arange(1, n_max+1), np.arange(1, n_max+1), psi_nm, shading='auto')
axs[1, 0].set_title('Transfer function (real part)')
axs[1, 0].pcolormesh(np.arange(1, n_max+1), np.arange(1, n_max+1), np.real(G), shading='auto')
axs[1, 1].set_title('Transfer function (imaginary part)')
axs[1, 1].pcolormesh(np.arange(1, n_max+1), np.arange(1, n_max+1), np.imag(G), shading='auto')
plt.show()
