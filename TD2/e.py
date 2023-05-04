import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parameters
a = 1  # Length of the membrane
N = 80  # Number of grid points in x and y direction
fps = 30  # Frames per second
t_end = 2  # End time of simulation
n_frames = int(t_end * fps)  # Number of frames in the animation

# Calculate the x and y coordinates of the grid points
x = np.linspace(0, a, N)
y = np.linspace(0, a, N)
X, Y = np.meshgrid(x, y)

# Calculate the eigenfrequencies and mode shapes
alpha_mn = np.zeros((N, N))
beta_mn = np.zeros((N, N))
omega_mn = np.zeros((N, N))
for m in range(1, N + 1):
    for n in range(1, N + 1):
        alpha_mn[m - 1, n - 1] = m * np.pi / a
        beta_mn[m - 1, n - 1] = n * np.pi / a
        omega_mn[m - 1, n - 1] = np.sqrt((alpha_mn[m - 1, n - 1] ** 2 + beta_mn[m - 1, n - 1] ** 2))

# Initialize the wave function
psi_nm = np.zeros((N, N))

# Create figure and axis objects
fig, ax = plt.subplots()

# Generate animation frames
frames = []
for i in range(n_frames):
    t = i / fps

    # Calculate the wave function at time t
    psi_nm = np.zeros((N, N))
    for m in range(1, N + 1):
        for n in range(1, N + 1):
            psi_nm += np.sin(alpha_mn[m - 1, n - 1] * X) * np.sin(beta_mn[m - 1, n - 1] * Y) * np.cos(omega_mn[m - 1, n - 1] * t)

    # Normalize the wave function and convert to uint8
    img = (psi_nm / np.max(psi_nm) * 255).astype(np.uint8)

    # Add the image to the list of frames
    frame = plt.imshow(img, cmap='gray', animated=True)
    frames.append([frame])

# Create the animation
ani = animation.ArtistAnimation(fig, frames, interval=1000/fps, blit=True)

# Save the animation as a GIF
ani.save('membrane.gif', writer='imagemagick')
