import numpy as np
import matplotlib.pyplot as plt

# Given parameters
a = 1  # m
b = 0.8  # m
mu_s = 0.262  # kg/m^2
T = 3990  # N/m

# Calculate the modal coefficients
N = 10  # number of modes
W_nm = np.zeros((N, N))
for n in range(1, N+1):
    for m in range(1, N+1):
        W_nm[n-1, m-1] = np.sqrt(4 * mu_s * T / (np.pi**4 * n**2 / a**2 + m**2 / b**2)**2)

# Visualize the modal coefficients
fig, ax = plt.subplots()
im = ax.imshow(W_nm, cmap='jet', origin='lower')
ax.set_title('Modal Coefficients $W_{nm}$')
ax.set_xlabel('m')
ax.set_ylabel('n')
cbar = ax.figure.colorbar(im)
cbar.ax.set_ylabel('$W_{nm}$', rotation=-90, va="bottom")
plt.show()

