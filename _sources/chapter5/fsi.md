# Fluid-Structure Interaction


Fluid-Structure Interaction (FSI) involves the coupled interaction between a fluid and a structure, where the motion of one influences the behavior of the other. The Finite Element Method (FEM) is a powerful tool for simulating FSI problems. This section provides an overview of the mathematical formulation for fluid-structure interaction using finite elements.

## Coupled Fluid-Structure Equations

The governing equations for fluid-structure interaction involve the Navier-Stokes equations for the fluid and the equations of motion for the structure. In a partitioned approach, the coupled system is given by:

$$
\begin{align}
    \text{Fluid Domain:} \quad \rho_f \left(\frac{\partial \mathbf{u}_f}{\partial t} + (\mathbf{u}_f \cdot \nabla)\mathbf{u}_f\right) &= -\nabla p_f + \mu_f \nabla^2 \mathbf{u}_f + \mathbf{f}_f + \mathbf{F}_s, \label{eq:fsi_fluid_momentum} \\
    \nabla \cdot \mathbf{u}_f &= 0, \label{eq:fsi_fluid_continuity} \\
    \text{Structure Domain:} \quad \rho_s \frac{\partial^2 \mathbf{u}_s}{\partial t^2} &= \nabla \cdot \mathbf{P}_s + \mathbf{f}_s, \label{eq:fsi_structure_motion}
\end{align}
$$

where:

- $\mathbf{u}_f$ is the fluid velocity,
- $p_f$ is the fluid pressure,
- $\rho_f$ is the fluid density,
- $\mu_f$ is the fluid dynamic viscosity,
- $\mathbf{f}_f$ is the fluid body force,
- $\mathbf{F}_s$ is the fluid-structure interaction force on the fluid by the structure,
- $\mathbf{u}_s$ is the structure displacement,
- $\rho_s$ is the structure density,
- $\mathbf{P}_s$ is the structure internal force,
- $\mathbf{f}_s$ is the structure external force.

## Coupling Conditions

For a fluid-structure interaction problem, coupling conditions at the fluid-structure interface are essential. These conditions ensure continuity of velocities, pressures, and forces between the fluid and structure. A common approach is to enforce kinematic and dynamic conditions:

$$
\begin{align}
    \text{Kinematic Condition:} \quad \mathbf{u}_f &= \mathbf{u}_s, \label{eq:fsi_kinematic_condition} \\
    \text{Dynamic Condition:} \quad \mathbf{P}_f \cdot \mathbf{n} &= \mathbf{P}_s \cdot \mathbf{n}, \label{eq:fsi_dynamic_condition}
\end{align}
$$

where $\mathbf{n}$ is the unit outward normal vector at the fluid-structure interface, and $\mathbf{P}_f$ and $\mathbf{P}_s$ are the fluid and structure stresses, respectively.

## Weak Formulation

The weak form of the coupled fluid-structure problem involves the variational formulation of the fluid and structure equations, incorporating the coupling conditions. For the fluid domain, the weak form is given by:

$$
\begin{equation}
    \int_{\Omega_f} \rho_f \left(\frac{\partial \mathbf{u}_f}{\partial t} + (\mathbf{u}_f \cdot \nabla)\mathbf{u}_f\right) \cdot \mathbf{v}_f \, d\Omega_f + \int_{\Omega_f} \mu_f \nabla \mathbf{u}_f : \nabla \mathbf{v}_f \, d\Omega_f = \int_{\Omega_f} (\mathbf{f}_f + \mathbf{F}_s) \cdot \mathbf{v}_f \, d\Omega_f,
\end{equation}
$$

where $\mathbf{v}_f$ is a test function belonging to the function space $H^1(\Omega_f)$, and $\Omega_f$ is the fluid domain.

For the structure domain, the weak form is given by:

$$
\begin{equation}
    \int_{\Omega_s} \rho_s \frac{\partial^2 \mathbf{u}_s}{\partial t^2} \cdot \mathbf{v}_s \, d\Omega_s + \int_{\Omega_s} \nabla \cdot \mathbf{P}_s \cdot \mathbf{v}_s \, d\Omega_s = \int_{\Omega_s} \mathbf{f}_s \cdot \mathbf{v}_s \, d\Omega_s,
\end{equation}
$$

where $\mathbf{v}_s$ is a test function belonging to the function space $H^1(\Omega_s)$, and $\Omega_s$ is the structure domain.

The coupling conditions, such as \eqref{eq:fsi_kinematic_condition} and \eqref{eq:fsi_dynamic_condition}, should be enforced during the assembly of the global matrices.

## References

The following references provide comprehensive coverage of the mathematical background for fluid-structure interaction using finite elements: {cite}`hron2006fluid`  {cite}`bathe2014finite` {cite}`quarteroni2017fluid`
