# the steady-state Navier-Stokes equations for incompressible fluids

The steady-state Navier Stokes problem for an incompressible fluid, with homogeneous Dirichlet boundary conditions (``no slip'' condition), is defined as

$$
\begin{equation}
\left\{
\begin{aligned}
        \left( \mathbf{u} \cdot \nabla \right) \mathbf{u} + \nabla p - 2 \nabla \cdot \left( \frac{1}{R_e} D\left( \mathbf{u} \right) \right) &= \mathbf{f} && \text{in $\Omega$},
  \\
        \nabla \cdot \mathbf{u} &= 0 && \text{in $\Omega$},
  \\
        \mathbf{u} &= 0 && \text{on $\partial\Omega$},
\end{aligned}
\right.
\end{equation}
$$

where $R_e$ is the Reynolds number and $D\left( \mathbf{u} \right) := \frac{1}{2}\left( \nabla \mathbf{u} + \nabla \mathbf{u}^T \right)$ is the strain rate tensor.
The weak formulation reads

$$
\begin{align}
  \text{find $\mathbf{u} \in W$ such that} \quad 
  a(\mathbf{u},\mathbf{v}) + b(\mathbf{u},\mathbf{v};\mathbf{u}) = l(\mathbf{v}) \quad
  \forall \mathbf{v} \in W
\end{align}
$$

where $W \subset \{ \mathbf{v} \in \mathbf{H}_0^1(\Omega): \nabla \cdot \mathbf{v} = 0 \}$ and

$$
\begin{align*}
  a(\mathbf{u}, \mathbf{v}) := \frac{2}{R_e} \int_{\Omega} D\left( \mathbf{u} \right) : D\left(\mathbf{v} \right) ~d\Omega,
  \quad
  b(\mathbf{u}, \mathbf{v}; \mathbf{w}) := \int_{\Omega} \left( \left( \mathbf{w} \cdot \nabla \right) \mathbf{u} \right) \cdot \mathbf{v} ~d\Omega,
  \quad
  l(\mathbf{v}) := \int_{\Omega} \mathbf{f} \cdot \mathbf{v} ~d\Omega.
\end{align*}
$$

