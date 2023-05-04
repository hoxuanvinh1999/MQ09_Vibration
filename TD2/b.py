import numpy as np
import matplotlib.pyplot as plt

# Dimensions de la membrane
a = 1.0   # Longueur
b = 0.8   # Largeur

# Propriétés de la membrane
mu_s = 0.262  # Masse surfacique
T = 3990.0    # Tension

# Densité de la membrane
rho = mu_s / (a*b)

# Coordonnées de l'excitation
x_e = a/5
y_e = b/6

# Fréquence d'excitation
f_e = 500.0   # Hz

# Calcul de la vitesse de propagation
c = np.sqrt(T/(rho * mu_s))

# Fonction de fréquence propre
def phi_n(x, n):
    return np.sin(n * np.pi * x / a)

def phi_m(y, m):
    return np.sin(m * np.pi * y / b)

# Calcul des fréquences propres et des coefficients modaux
n_max = 10
m_max = 10

k_nm = np.zeros((n_max, m_max))
A_nm = np.zeros((n_max, m_max))

for n in range(1, n_max+1):
    for m in range(1, m_max+1):
        k_nm[n-1, m-1] = np.pi * np.sqrt((n/a)**2 + (m/b)**2)
        A_nm[n-1, m-1] = 1/(rho * k_nm[n-1, m-1]**2) * np.sin(n * np.pi * x_e / a) * np.sin(m * np.pi * y_e / b)

# Calcul de la fonction de transfert pour chaque mode
n_values = np.arange(1, n_max+1)
m_values = np.arange(1, m_max+1)
H_nm = np.zeros((n_max, m_max))

omega = 2 * np.pi * f_e

for n in n_values:
    for m in m_values:
        H_nm[n-1, m-1] = A_nm[n-1, m-1] / (rho * k_nm[n-1, m-1] * phi_n(x_e, n) * phi_m(y_e, m) * omega**2 - T)

# Affichage de la fonction de transfert sous forme de carte de chaleur
plt.imshow(H_nm, cmap='plasma')
plt.colorbar()
plt.title('Fonction de transfert pour une excitation à {} Hz en ({}, {})'.format(f_e, x_e, y_e))
plt.xlabel('n')
plt.ylabel('m')
plt.show()