## Overview

The model solves:
```math
\frac{\partial N}{\partial t} + \frac{\partial (C_x N)}{\partial x} + \frac{\partial (C_y N)}{\partial y} + \frac{\partial (C_\sigma N)}{\partial \sigma} + \frac{\partial (C_\theta N)}{\partial \theta} = \frac{S}{\sigma}
```

Where:
- $t$ : time
- $x, y$ : spatial coordinate
- $\sigma$ : intrinsic wave frequency
- $\theta$ : wave direction
- $C_i$ : propagation velocities
- $N(x, y, \sigma, \theta, t)$ : wave action density
- $S$ : net source and sink term

This equation is implemented in a modular 5D solver inspired by spectral wave model [WAVEWATCH III](https://github.com/NOAA-EMC/WW3)

## Project Structure
```markdown
wavemod/
├── src/
|   ├── grid.f90        <- Grid initialization for x, y, sigma, theta
|   ├── windmod.f90     <- Wind initialization
|   ├── physics.f90     <- Wind input (S_in) and dissipation (S_ds) (RHS)
|   ├── waveprop.f90    <- wave propagation (spatial and intra-spectral) (2nd, 3rd, 4th, and 5th terms in LHS)
|   ├── io.f90          <- Input/output utilities
|   ├── params.f90      <- Wrapper for params calculation
|   ├── main.f90        <- Time integration and main loop (1st term in LHS)
|   └── Makefile        <- Build configuration
├── py/
|   ├── input.ipynb     <- Wind generation
|   └── output.ipynb    <- Output visualization
├── input/*.txt         <- Input fields (wind)
└── output/             <- Output fileds (Hs and Tm)
    ├── diag/*.txt      <- Output Hs and Tm
    └── spc/*.txt       <- Output in spectral form (1 grid)
```

## Requirements

- Fortran 90+ compiler (for now only tested using `gfortran`)
- `make` utility
- Python (for pre and post processing)

## Features

- 2D spatial grid $(x, y)$
- 1D spectral bins:
    - Frequency $(\sigma)$: logaritmically spaced
    - Direction $(\theta)$: uniformly spaced
- Explicit time stepping (Euler forward)
- First-order upwind advection (spatial and intra wave propagation)
- Modular source terms
- ASCII-based I/O (easy to inspect and debug)

## Tunable constants and model configurations

The following constants are defined in `src/physics.f90` may be modified:

| Constant      | Description                   | Default   |
| -             | -                             | -         |
| `A_in`        | Wind input growth factor      | 1.90      |
| `B_ds`        | Dissipation coefficient       | 0.0002    |
| `m_ds`        | Dissipation exponent          | 4.0       |
| `sigma_peak`  | Peak freq for whitecapping    | 0.20      |
| `g`           | Gravitational acceleration    | 9.81      |
|               |                               |           |

The following model configurations are defined in `src/grid.f90` may be modified:
| Config param  | Description                   | Default   |
| -             | -                             | -         |
| `nx`, `ny`    | Number of spatial grids       | 30 x 30   |
| `dx`, `dy`    | Spatial resolution            | 1.0       |
| `dt`          | Time step (in minutes)        | 10.0      |
| `tmax`        | Total simulation time         | 3600 * 1  |
| `nt`          | Number of time steps          | `tmax/dt` |
|               |                               |           |
| `nsigma`      | Number of freq bins           | 30        |
| `sigma_min`   | Minimum freq                  | 0.04      |
| `sigma_max`   | Maximum freq                  | 1         |
| `sigma(i)`    | Freq bins (log-distributed)   | $\sigma_i = \sigma_{min} \cdot \exp(d_{log(\sigma)} \cdot(i-1))$
| `dsigma(i)`   | Bin width for frequency (computed from $\sigma$) | Forward difference (last repeated) |
|               |                               |           |
| `ntheta`      | Number of directional bins    | 36        |
| `theta(i)`    | Direction angle (radians)     | Uniform ($0$ to $2\pi$) |
|               |                               |           |

## Outputs
- `Hs`: signficant wave height
- `Tm`: mean wave period
- Outputs are stored in `output` directory

#### Wave Spectra comparison
The wave spectra for a specific grid point is compared in `py/output.py` file with Pierson-Moskowitz and JONSWAP spectrum.

