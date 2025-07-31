# Linear Elasticity

## The PDE problem
The governing equations for small elastic deformations of a body $\Omega$ can be expressed as:

$$
\begin{aligned}
    -\nabla \cdot \sigma(u) &= f & \text{in } \Omega \\
    \sigma(u) &= C : \epsilon(u) & \text{in } \Omega \\
\end{aligned}
$$

where : 
 - $\Omega$ is the domain of interest,
 - $u$ is the displacement vector field,
 - $\sigma(u)$ is a second-order symmetric tensor, called stress tensor, 
 - $f$ represents the load vector field (force per unit volume), 
 <!-- $\kappa$ and $\mu$ are Lamé's elasticity parameters for the material, $I$ denotes the identity tensor,  -->
 - $\epsilon(u)$ is a second-order symmetric tensor, called strain tensor and by definition $\epsilon(u) := \frac{1}{2}(\nabla u + (\nabla u)^T)$.
 - $C$ is the elasticity tensor, which relates stress and strain. 

Then, the strong formulation of linear elasticity is :
Find $u \in V$ such that

$$
\begin{aligned}
    -\nabla \cdot \sigma (u) &= f & \text{in } & \Omega \\
    u &= 0 & \text{on } & \partial \Omega_D \\
    \sigma(u) \cdot n &= g_T & \text{on } & \partial \Omega_T
\end{aligned}
$$

where : 
- $V = \{ v \in (H^1(\Omega))^3 : v = 0 \text{ on } \partial \Omega_D \}$,

- $g_T$ is the traction vector on the part $\partial \Omega_T$ of the boundary,

- $n$ is the outward normal vector on the boundary $\partial \Omega$,

- $\partial \Omega_D \cap \partial \Omega_T = \emptyset$ and $\partial \Omega = \partial \Omega_D \cup \partial \Omega_T$. 
- $\partial \Omega_D$ is the Dirichlet boundary condition where the displacement is fixed to zero, and $\partial \Omega_T$ is the Neumann boundary condition where the traction is prescribed.

## The Variational Formulation
The variational formulation of the linear elasticity equations involves forming the inner product of the PDE with a vector test function $v \in V$ and integrating over the domain $\Omega$. This yields:

$$
\int_{\Omega} - \nabla \cdot \sigma(u) \cdot v \, \mathrm{d} x = \int_{\Omega} f \cdot v \, \mathrm{d} x
$$

Integrating the term $\nabla \cdot \sigma(u) \cdot v$ by parts, considering boundary conditions, we obtain:

$$
\int_{\Omega} \sigma(u) : \nabla v \, \mathrm{d} x = \int_{\Omega} f \cdot v \, \mathrm{d} x + \int_{\partial \Omega_T} g_T \cdot v \, \mathrm{d} s
$$

By using the symmetry of the stress tensor $\sigma$ and its definition from $(2)$, we can notice that :

$$
\begin{aligned}
\int_{\Omega} \sigma(u) : \nabla v \, \mathrm{d} x &= \int_{\Omega} \sigma(u) : \epsilon(v) \, \mathrm{d} x = \int_{\Omega} C : \epsilon(u) : \epsilon(v) \, \mathrm{d} x \\ &= \int_{\Omega} \epsilon(u) : C : \epsilon(v) \, \mathrm{d} x
\end{aligned}
$$

This leads to the following variational formulation:

$$
\boxed{
\begin{aligned}
&\text{Find } u \in V \text{ such that:} \\
&\qquad a(u, v) = L(v) \quad \forall v \in V
\end{aligned}
}
$$

with

$$
\begin{aligned}
&a : 
\begin{cases}
V \times V \rightarrow \mathbb{R} \\
(u, v) \longmapsto \int_{\Omega} \epsilon(u) : C : \epsilon(v) \, \mathrm{d} x
\end{cases} \\[0.3cm]
&L : 
\begin{cases}
V \rightarrow \mathbb{R} \\
v \longmapsto \int_{\Omega} f \cdot v \, \mathrm{d} x + \int_{\partial \Omega_T} g_T \cdot v \, \mathrm{d} s
\end{cases}
\end{aligned}
$$

## Isotropic Materials
For isotropic materials, the elasticity tensor $C$ can be expressed in terms of the Lamé parameters $\lambda$ and $\mu$ as follows:

$$
C := \lambda (\nabla \cdot u) I_3 + 2\mu \epsilon(u)
$$
Then, the stress tensor can be expressed as:
$$\sigma(u) = \lambda (\nabla \cdot u) I_3 + 2\mu \epsilon(u)$$


This leads to the variational formulation:

$$
\boxed{
\begin{aligned}
&\text{Find } u \in V \text{ such that:}
\\
&\qquad a(u, v) = L(v) \quad \forall v \in V
\end{aligned}
}
$$

with

$$
\begin{aligned}
&a :\begin{cases}
V \times V \rightarrow \mathbb{R} \\
(u, v) \longmapsto \int_{\Omega} \sigma(u) : \epsilon(v) \, \mathrm{d} x
\end{cases} \\[0.3cm]
&L :\begin{cases}
V \rightarrow \mathbb{R} \\
v \longmapsto \int_{\Omega} f \cdot v \, \mathrm{d} x + \int_{\partial \Omega_T} g_T \cdot v \, \mathrm{d} s
\end{cases}
\end{aligned}
$$

With this formulation, the problem is well-posed under the assumption that the material is isotropic and the boundary conditions are properly defined. While $\frac{\lambda}{\mu}$ is not too large (typically $\frac{\lambda}{\mu} \leq 10^4$), the problem remains well-posed numerically. However, as $\frac{\lambda}{\mu}$ increases, the problem can become ill-posed, leading to numerical difficulties in finding a solution. The first notebook of this chapter illustrates the case of isotropic materials with $\frac{\lambda}{\mu} \leq 10^4$ and the second notebook is trying to illustrate the case of isotropic materials with $\frac{\lambda}{\mu} > 10^4$.

<p align="center">

| Material  | $\lambda$ (GPa) | $\mu$ (GPa) |
|-----------|:--------------:|:-----------:|
| Steel     | 120            | 80          |
| Concrete  | 17             | 14          |
| Rubber    | 0.16           | 0.00033     |

</p>
