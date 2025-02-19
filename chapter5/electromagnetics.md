# Electromagnetics


Electromagnetic problems are commonly described by Maxwell's equations, which govern the behavior of electric and magnetic fields. The Finite Element Method (FEM) provides a powerful numerical approach for solving these equations in complex geometries. This section provides a concise overview of the mathematical formulation for electromagnetic problems using finite elements.

## Maxwell's Equations

The time-harmonic Maxwell's equations in free space are given by:

$$
\begin{align}
    \nabla \times \mathbf{E} &= -j\omega \mu \mathbf{H}, \label{eq:maxwell1} \\
    \nabla \times \mathbf{H} &= j\omega \varepsilon \mathbf{E}, \label{eq:maxwell2} \\
    \nabla \cdot \mathbf{B} &= 0, \label{eq:maxwell3} \\
    \nabla \cdot \mathbf{D} &= 0, \label{eq:maxwell4}
\end{align}
$$

where $\mathbf{E}$ is the electric field, $\mathbf{H}$ is the magnetic field, $\mathbf{B}$ is the magnetic flux density, $\mathbf{D}$ is the electric displacement field, $\omega$ is the angular frequency, $\mu$ is the permeability, and $\varepsilon$ is the permittivity.

## Weak Formulation

The weak form of Maxwell's equations is obtained by multiplying each equation with suitable test functions and integrating over the domain $\Omega$. For the electric field $\mathbf{E}$, the weak form is given by:

$$
\begin{equation}
    \int_{\Omega} \nabla \times \mathbf{E} \cdot \nabla \times \mathbf{v} \, d\Omega + \omega^2 \mu \varepsilon \mathbf{E} \cdot \mathbf{v} \, d\Omega = 0,
\end{equation}
$$

where $\mathbf{v}$ is a test function belonging to the function space $H(\text{curl}; \Omega)$.

Similarly, for the magnetic field $\mathbf{H}$, the weak form is given by:

$$
\begin{equation}
    \int_{\Omega} \nabla \times \mathbf{H} \cdot \nabla \times \mathbf{u} \, d\Omega - \omega^2 \mu \varepsilon \mathbf{H} \cdot \mathbf{u} \, d\Omega = 0,
\end{equation}
$$

where $\mathbf{u}$ is a test function belonging to the function space $H(\text{curl}; \Omega)$.

## Finite Element Discretization

To apply the finite element method, the domain $\Omega$ is discretized into elements, and the electromagnetic fields are approximated using piecewise basis functions. The discretized weak form leads to a system of linear equations, which can be solved numerically to obtain the finite element solution.

## References

The following references provide in-depth coverage of the mathematical background for electromagnetism using finite elements: {cite}`jin2015finite`  {cite}`bossavit1998whitney` {cite}`monk2003finite`

