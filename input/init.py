import numpy as np

# Redefine constants after kernel reset
nsteps = 72
nx, ny = 10, 10
nsigma, ntheta = 30, 36

# Generate realistic wind field: sinusoidal variation + small noise
u10 = np.zeros((nsteps, ny, nx))
for n in range(nsteps):
    wind_base = 5 + 5 * np.sin(2 * np.pi * n / nsteps)
    u10[n] = wind_base + 0.5 * np.random.randn(ny, nx)

# Generate boundary spectrum N_bdy: time-varying flat low energy
N_bdy = np.zeros((nsteps, nx, ny, nsigma, ntheta))
for n in range(nsteps):
    for k in range(nsigma):
        for l in range(ntheta):
            N_bdy[n, :, :, k, l] = 0.01 + 0.005 * np.sin(2 * np.pi * n / nsteps)

# Save wind forcing
wind_file = "/home/yothunder/fort/wvmod/input/u10.txt"
with open(wind_file, "w") as f:
    for n in range(nsteps):
        for j in range(ny):
            for i in range(nx):
                f.write(f"{u10[n, j, i]:.4f}\n")

# Save boundary spectrum
boundary_file = "/home/yothunder/fort/wvmod/input/bound.txt"
with open(boundary_file, "w") as f:
    for n in range(nsteps):
        for l in range(ntheta):
            for k in range(nsigma):
                for j in range(ny):
                    for i in range(nx):
                        f.write(f"{N_bdy[n, i, j, k, l]:.5f}\n")

wind_file, boundary_file
