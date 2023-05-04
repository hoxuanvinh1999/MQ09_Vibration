import numpy as np
import matplotlib.pyplot as plt

# Constants
a = 1.0  # length of membrane (m)
b = 0.8  # width of membrane (m)
mu_s = 0.262  # mass per unit area (kg/m^2)
T = 3990  # tension (N/m)

# Frequencies to evaluate
freqs = np.linspace(1, 500, 5000)  # Hz

# Calculate wavenumbers
k_x = np.pi / a * np.arange(1, 21)
k_y = np.pi / b * np.arange(1, 17)

# Calculate modal mass
M_mn = mu_s * a * b / np.outer(k_x ** 2, k_y ** 2)

# Calculate natural frequencies
omega_mn = np.outer(k_x, np.ones_like(k_y)) ** 2 + np.outer(np.ones_like(k_x), k_y) ** 2
omega_mn = np.sqrt(T / mu_s * omega_mn)
f_mn = omega_mn / (2 * np.pi)  # convert to Hz

# Calculate modal response to uniform excitation
F_mn = 4 * np.pi ** 4 * mu_s * a * b * (freqs[:, None, None] ** 2 - omega_mn ** 2) ** -2
psi_nm = F_mn / M_mn

# Calculate transfer function
H = np.sum(psi_nm, axis=(1, 2)) / np.sum(F_mn, axis=(1, 2))

# Plot results
plt.figure()
plt.plot(freqs, np.abs(H))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Frequency response of rectangular membrane')
plt.show()
