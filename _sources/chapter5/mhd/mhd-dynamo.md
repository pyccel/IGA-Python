# MHD Dynamo Theory

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

3. MHD Induction Equation:
   $$
   \frac{\partial \mathbf{B}}{\partial t} = \nabla \times (\mathbf{v} \times \mathbf{B} - \eta \nabla \times \mathbf{B})
   $$
   where $\eta$ is the magnetic diffusivity.

4. MHD Ohm's Law:
   $$
   \mathbf{J} = \sigma (\mathbf{E} + \mathbf{v} \times \mathbf{B})
   $$
   with $\sigma$ as the electrical conductivity and $\mathbf{E}$ as the electric field.

5. MHD Dynamo Equation:
   $$
   \frac{\partial \mathbf{B}}{\partial t} = \nabla \times (\mathbf{v} \times \mathbf{B} - \eta \nabla \times \mathbf{B} + \alpha \mathbf{B})
   $$
   introducing the dynamo term ($\alpha \mathbf{B}$) representing the generation of magnetic field.

