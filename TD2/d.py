import numpy as np
import matplotlib.pyplot as plt

# Constants
a = 1.0  # m
b = 0.8  # m
mu_s = 0.262  # kg/m^2
T = 3990  # N/m
rho = 7850  # kg/m^3

# Modes
n_max = 5
m_max = 5

# Frequencies
f_e = 500  # Hz
omega = 2 * np.pi * f_e

# Excitation position
x_e = a / 5
y_e = b / 6

# Grid
N = 100
M = 80
X, Y = np.meshgrid(np.linspace(0, a, N), np.linspace(0, b, M))

# Coefficients
k_nm = np.zeros((n_max, m_max))
A_nm = np.zeros((n_max, m_max))
W_nm = np.zeros((n_max, m_max))

for n in range(1, n_max + 1):
    for m in range(1, m_max + 1):
        k_nm[n-1, m-1] = np.sqrt((n*np.pi/a)**2 + (m*np.pi/b)**2)
        A_nm[n-1, m-1] = (np.sin(n * np.pi * x_e / a)
                          * np.sin(m * np.pi * y_e / b))
        W_nm[n-1, m-1] = np.sqrt(mu_s / (A_nm[n-1, m-1] * T)) * k_nm[n-1, m-1] / np.sqrt(rho)

# Deformed shape
u = np.zeros_like(X)

for n in range(1, n_max + 1):
    for m in range(1, m_max + 1):
        u += W_nm[n-1, m-1] * np.sin(n * np.pi * X / a) * np.sin(m * np.pi * Y / b)

# Plotting
fig, ax = plt.subplots()
c = ax.contourf(X, Y, u, cmap='coolwarm')
plt.colorbar(c)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('Deformée opérationnelle')
plt.show()
