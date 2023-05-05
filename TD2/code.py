# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:47:38 2023

@author: GLUTTONY
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 10:51:28 2019
"""
from pylab import *
from scipy.io.wavfile import write
import winsound
import matplotlib.pyplot as plt
import matplotlib.animation as animation

## Données géométriques
a = 1.0
b = 0.8

## Propriétés de la membrane
T = 3990.0        # Tension membrane
mu = 0.262        # masse surfacique

c = sqrt(T/mu)

## Decomposition modale
N_modes = 200
N = 19
n = arange(1, N+1)
n = reshape(n, (N, 1))
M = 21
m = arange(1, M+1)
m = reshape(m, (M, 1))

omeganm = c*sqrt((n*pi/a)**2 + (m.T*pi/b)**2)

omega = sort(ravel(omeganm))
idx = argsort(ravel(omeganm))

f = 1/(2*pi)*omega

print('10 premières fréquences (Hz)')
print(f[:10])

Nx = 200
Ny = 200

x = linspace(0, a, Nx)
y = linspace(0, b, Ny)
[X, Y] = meshgrid(x, y)

Phinm = zeros((N, M, Nx, Ny))
for i in range(N):
    for j in range(M):
        Phinm[i, j, :, :] = sin(n[i]*pi*X/a)*sin(m[j]*pi*Y/b)

Phi = reshape(Phinm, (N*M, Nx, Ny))

Phi = Phi[idx, :, :]

## Affichage des déformées modales
figure()
for i in range(4):
    subplot(2, 2, i+1)
    contourf(X, Y, Phi[i, :, :])
    xlabel('x (m)')
    ylabel('y (m)')
    title('Déformée modale ' + str(i+1))

show()


## REPONSE FORCEE
fe = 500  #Hz
we = 2*pi*fe
xe = a/5
ye = b/6
ksi = 0.001

Mnm = mu*a*b/4

## a) Force modale

Fnm = zeros((N, M))


for i in range(N):
    for j in range(M):
        Fnm[i,j] = 4*Mnm*ksi*we**2*Phinm[i,j,int(xe*Nx/a),int(ye*Ny/b)]

figure()
imshow(abs(Fnm),   extent = [0.5, N+0.5, 0.5, M+0.5], origin = 'lower')
colorbar()
title('Force modale')
xlabel('n')
ylabel('m')


## b) Fonction de transfert
Hnm = zeros((N, M), dtype=complex)

for n in range(1, N+1):
    for m in range(1, M+1):
        Hnm[n-1,m-1] = (1/(mu*a*b*we**2 - pi**4*T/(4*a**2*b**2))) * (4/(pi**2*(n**2/a**2 + m**2/b**2))) 

figure()
imshow(abs(Hnm),   extent = [0.5, N+0.5, 0.5, M+0.5], origin = 'lower')
xlabel('n')
ylabel('m')
title('Fonction de transfert')

## c) Coefficients modaux
Wnm = Hnm*Fnm

figure()
imshow(abs(Wnm),   extent = [0.5, N+0.5, 0.5, M+0.5], origin = 'lower')
xlabel('n')
ylabel('m')
title('Coefficients modaux')

show()

## d) Déformée opérationnelle
u = zeros((Nx, Ny), dtype=complex)
for n in range(1, N+1):
    for m in range(1, M+1):
        u += real(Wnm[n-1, m-1])*Phinm[n-1, m-1, :, :]

u = real(u)

figure()
contourf(X, Y, u)
colorbar()
xlabel('x (m)')
ylabel('y (m)')
title('Déformée opérationnelle')

show()

## e) 
from matplotlib.animation import FuncAnimation

# Paramètres de l'animation
duration = 5  # durée de l'animation en secondes
fps = 30  # nombre d'images par seconde

# Fonction qui calcule la déformée à chaque instant
def animate(t):
    # Décalage de la position d'excitation au cours du temps
    xe = a/5 + 0.05*a*cos(2*pi*fe*t)
    ye = b/6 + 0.05*b*sin(2*pi*fe*t)
    
    # Calcul de la déformée opérationnelle
    u = zeros((Nx, Ny), dtype=complex)
    for n in range(1, N+1):
        for m in range(1, M+1):
            # Calcul de la force modale
            Fnm = 4*Mnm*ksi*we**2*sin(n*pi*xe/a)*sin(m*pi*ye/b)
            
            # Calcul de la fonction de transfert
            Hnm = (1/(mu*a*b*we**2 - pi**4*T/(4*a**2*b**2))) * (4/(pi**2*(n**2/a**2 + m**2/b**2)))
            
            # Calcul du coefficient modal
            Wnm = Hnm*Fnm
            
            # Calcul de la déformée modale
            Phinm = sin(n*pi*X/a)*sin(m*pi*Y/b)
            
            # Calcul de la déformée opérationnelle
            u += real(Wnm)*Phinm
    
    # Affichage de la déformée opérationnelle
    ax.clear()
    ax.contourf(X, Y, u, levels=50)
    ax.set_title('Déformée opérationnelle\nt = {:.2f} s'.format(t))
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_xlim([0, a])
    ax.set_ylim([0, b])
    
    return ax,

# Création de la figure
fig, ax = plt.subplots()

# Animation de la déformée opérationnelle
animation = FuncAnimation(fig, animate, frames=np.linspace(0, duration, duration*fps), blit=True)

# Sauvegarde de l'animation en format gif
animation.save('deformee_operationnelle.gif', writer='imagemagick', fps=fps)

# Affichage de l'animation
plt.show()














