# Magnetic Drug Targeting

Magnetic Drug Targeting involves the use of magnetic fields to guide drug-carrying particles to a specific target within the body. This application combines principles of magnetohydrodynamics (MHD) with drug delivery. Below are the mathematical models and example test cases for Magnetic Drug Targeting:

## Mathematical Models

1. MHD Continuity Equation:
   $$
   \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0
   $$
   where $\rho$ is the fluid density, and $\mathbf{v}$ is the fluid velocity.

2. MHD Momentum Equation:
   $$
   \rho \frac{D\mathbf{v}}{Dt} = -\nabla p + \rho \mathbf{g} + \nabla \cdot \boldsymbol{\tau} + \mathbf{J} \times \mathbf{B} + \nu \nabla^2 \mathbf{v}
   $$
   incorporating the viscous term ($\nu \nabla^2 \mathbf{v}$) for fluid viscosity.

3. Magnetic Particle Motion Equation:
   $$
   m_p \frac{d\mathbf{u}_p}{dt} = \mathbf{F}_p + \mathbf{F}_{\text{magnetic}}
   $$
   where $m_p$ is the particle mass, $\mathbf{u}_p$ is the particle velocity, $\mathbf{F}_p$ is the sum of external forces, and $\mathbf{F}_{\text{magnetic}}$ is the magnetic force.

4. Magnetic Force Equation:
   $$
   \mathbf{F}_{\text{magnetic}} = \nabla (\mathbf{m} \cdot \mathbf{B})
   $$
   where $\mathbf{m}$ is the magnetic moment of the drug-carrying particle.

