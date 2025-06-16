## Overview

The model solves:
$$
\frac{\partial N}{\partial t} + \frac{\partial (C_x N)}{\partial x} + \frac{\partial (C_y N)}{\partial y} + \frac{\partial (C_\sigma N)}{\partial \sigma} + \frac{\partial (C_\theta N)}{\partial \theta} = \frac{S}{\sigma}
$$

Where:
- $(N(x, y, \sigma, \theta, t))$ is the wave action density
- $(C_i)$ are the propagation velocities
- $(S)$ represents source/sink terms (wind input & dissipation)

## Project Structure
```markdown
wavemod/
├── grid.f90        <- Grid initialization for x, y, sigma, theta
├── windmod.f90     <- Wind initialization
├── physics.f90     <- Wind input (S_in) and dissipation (S_ds) (RHS)
├── waveprop.f90    <- wave propagation (spatial and intra-spectral) (2nd, 3rd, 4th, and 5th terms in LHS)
├── io.f90          <- Input/output utilities
├── params.f90      <- Wrapper for params calculation
├── main.f90        <- Time integration and main loop (1st term in LHS)
├── Makefile        <- Build configuration
├── input/          <- Input fields (wind)
└── output/         <- Output fileds (Hs and Tm)
```

## Requirements

- Fortran 90+ compiler (for now only tested using `gfortran`)
- Make utility
- Python (for pre and post processing)

## Features

- 2D spatial grid $(x, y)$
- 1D logarithmic frequency bins $(\sigma)$
- 1D uniform direction bins $(\theta)$
- Explicit time stepping
- First-order upwind advection (spatial and intra wave propagation)
- Modular source term implementation
- ASCII-based I/O