# Linear Elasticity

## The PDE problem
The governing equations for small elastic deformations of a body $\Omega$ can be expressed as:

$$
\begin{align}
    -\nabla \cdot \sigma(u) &= f & \text{in } \Omega \\
    \sigma(u) &= \kappa \text{tr}(\epsilon(u))I + 2 \mu \epsilon(u)
\end{align}
$$

where $\sigma$ is the stress tensor, $f$ represents the body force per unit volume, $\kappa$ and $\mu$ are Lamé's elasticity parameters for the material, $I$ denotes the identity tensor, $\epsilon$ is the symmetric strain tensor and the displacement vector field is denoted by $u$. $\epsilon := \frac{1}{2}(\nabla u + (\nabla u)^T)$



By substituting $\epsilon(u)$ into $\sigma$, we obtain:

$$
\sigma(u) = \kappa (\nabla \cdot u)I + \mu(\nabla u + (\nabla u)^T)
$$

Then, the strong formulation of linear elasticity is :
Find $u \in V$ such that

$$
\begin{align}
    -\nabla \cdot \sigma (u) &= f & \text{in } & \Omega \\
    u &= 0 & \text{on } & \partial \Omega_D \\
    \sigma(u) \cdot n &= g_T & \text{on } & \partial \Omega_T
\end{align}
$$

where : 
- $V = \{ v \in (H^1(\Omega))^3 : v = 0 \text{ on } \partial \Omega_D \}$,

- $g_D$ is the Dirichlet boundary condition on the part $\partial \Omega_D$ of the boundary,

- $g_T$ is the traction vector on the part $\partial \Omega_T$ of the boundary,

- $n$ is the outward normal vector on the boundary $\partial \Omega$,

- $\partial \Omega_D \cap \partial \Omega_T = \emptyset$ and $\partial \Omega = \partial \Omega_D \cup \partial \Omega_T$.
