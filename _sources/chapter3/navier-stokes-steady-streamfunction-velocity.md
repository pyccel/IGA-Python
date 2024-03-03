# The streamfunction-velocity formulation of the steady-state Navier-Stokes equations for incompressible fluids
*Author: Ahmed Ratnani*

When $\Omega$ is a simply connected 2D domain, there exists a unique function $\psi$ such that $\mathbf{u} = \boldsymbol{\nabla} \times \psi:= \left( \partial_y \psi, - \partial_x \psi \right)$;
substituting this expression for $\mathbf{u}$ into \eqref{eq:steady-navier-stokes} leads to the so-called ``streamfunction-velocity formulation'' of the steady-state Navier-Stokes equations for an incompressible fluid.
For more details we refer the reader to \cite{TAGLIABUE2014277}.
The variational formulation is then given by:

$$
\begin{align}
  \text{find $\psi \in V$ such that} \quad 
  a(\psi,\phi) + b(\psi,\phi;\psi) = l(\phi)
  \quad \forall \phi \in V,
\end{align}
$$

where $V \subset H^2(\Omega)$ and  

$$
\begin{align*}
  a(\psi, \phi) := \int_{\Omega} D\left( \boldsymbol{\nabla} \times \psi \right) : D\left(\boldsymbol{\nabla} \times \phi \right) ~d\Omega,  
  \quad
  b(\psi, \phi; \xi) := \int_{\Omega} \left( \left(\boldsymbol{\nabla} \times \xi \cdot \nabla \right) \boldsymbol{\nabla} \times \psi \right) \cdot \boldsymbol{\nabla} \times \phi ~d\Omega,
  \quad
  l(\phi) := \int_{\Omega} \mathbf{f} \cdot \boldsymbol{\nabla} \times \phi ~d\Omega.
\end{align*}
$$

%As in the previous example, one can use a Picard method, by considering $\xi$ as a \texttt{Field} in Eq. \ref{eq:steady-streamfunction-velocity-vf-a}, which will lead to the Python code described in \ref{code:steady-streamfunction-velocity-picard}. %, or consider $a$ as a \texttt{LinearForm} where $\xi = \psi$ are \texttt{Field} objects, and use the \texttt{NewtonIteration} constructor to get a \texttt{BilinearForm}, as described in the Python code \ref{code:steady-streamfunction-velocity-newton}.


## Linearization of the streamfunction-velocity linear form

In the sequel, we use SymPDE to linearize the streamfunction-velocity formulation, starting from a \texttt{LinearForm}.
As presented in the previous example, we shall write our weak formulation as a linear form where the velocity $\mathbf{u}$ (\textit{i.e.} \texttt{u}) is considered as a \textit{field}, which is implemented in Python code \ref{code:steady-navier-stokes-nl}.
Using the function \texttt{linearize} allows us to construct the bilinear form associated to the linearization of the given linear form, around the field \texttt{psi}, see Python code \ref{code:steady-streamfunction-velocity-lin}.
To get an \texttt{Equation}, the user can add the left-hand-side for the Navier-Stokes equation to the linear form, the appropriate essential boundary conditions, then use the \texttt{NewtonIteration} as presented in the previous example.
