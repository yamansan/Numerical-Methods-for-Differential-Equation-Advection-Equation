# Numerical Simulation of Advection and Burgers' Equations

## Overview
This repository provides an in-depth numerical simulation of the **Advection Equation** and the **Burgers' Equation**, two fundamental models in fluid dynamics and traffic flow. These simulations employ advanced finite difference methods to ensure accuracy and stability while preserving the physical properties of the system, such as volume and mass conservation.

Our approach utilizes:

- **Upwind Scheme**
- **Forward Difference Scheme**
- **Central Difference Scheme**
- **Min-Mod Limiter**

By leveraging these techniques, the project achieves robust and accurate numerical solutions suitable for modeling real-world phenomena.

## Mathematical Background

### Advection Equation
The advection equation describes the transport of a scalar quantity (e.g., mass or energy) in a flow field. The one-dimensional form is given by:

$$
\frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0,
$$

where:
- \( u(x,t) \) is the transported quantity,
- \( c \) is the constant advection velocity.

### Burgers' Equation
The Burgers' equation models viscous flow and traffic dynamics. The one-dimensional form is:

$$
\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} = \nu \frac{\partial^2 u}{\partial x^2},
$$

where:
- \( u(x,t) \) is the velocity field,
- \( \nu \) is the kinematic viscosity.

## Numerical Methods

### Finite Difference Schemes
The project implements the following finite difference methods:

1. **Upwind Scheme**:
   - Introduces numerical diffusion to stabilize the solution.
   - Suitable for hyperbolic problems with shock waves.
   - Discretized as:
     $$
     \frac{u^{n+1}_i - u^n_i}{\Delta t} + c \frac{u^n_i - u^n_{i-1}}{\Delta x} = 0.
     $$

2. **Forward Difference Scheme**:
   - Uses forward time differencing for time derivatives:
     $$
     \frac{u^{n+1}_i - u^n_i}{\Delta t} = -c \frac{u^n_{i+1} - u^n_{i-1}}{2\Delta x}.
     $$

3. **Central Difference Scheme**:
   - Provides second-order accuracy in space:
     $$
     \frac{\partial u}{\partial x} \approx \frac{u_{i+1} - u_{i-1}}{2\Delta x}.
     $$

### Min-Mod Limiter
To address the challenges of shock capturing and oscillations, the **Min-Mod Limiter** is employed. It ensures the preservation of mass and volume by preventing spurious oscillations:

$$
\phi(r) = \text{minmod}(r, 1) = \begin{cases}
    r & \text{if } 0 \leq r \leq 1, \\
    1 & \text{if } r > 1, \\
    0 & \text{if } r < 0.
\end{cases}
$$

This limiter adapts the slope in regions with sharp gradients, balancing accuracy and stability.

## Applications

### Fluid Dynamics
The simulation captures the dynamics of wave propagation and shock formation in fluids. It preserves the mass and volume of the flow, making it suitable for modeling compressible and incompressible flows.

### Traffic Flow
The Burgers' equation is applied to simulate traffic dynamics, including congestion and wave-like propagation of density.

## Getting Started

### Prerequisites
- Python 3.8+
- NumPy
- Matplotlib



