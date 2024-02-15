# The Poisson problem using Nitsche method

In this section, we consider a domain $\Omega$ which a union of multiple subdomain, *i.e.*

$$
\Omega := \bigcup_{i=1}^n \Omega_i, \quad \Omega_i \bigcap \Omega_j = \emptyset, ~ i \neq j  
$$

We shall denote by $\mathcal{I}$ the set of all internal interfaces of $\Omega$.

We shall also need the following operators

- The jump of the function $u$, defined as $[\![ u ]\!] := u|\_{\Omega_{i_1}} -  u|\_{\Omega_{i_2}}$ for two adjacent subdomains $\Omega_{i_1}$ and $\Omega_{i_2}$
- The average of the function $u$, defined as $\{u\} := \frac{1}{2} \left( u|\_{\Omega_{i_1}} +  u|\_{\Omega_{i_2}} \right)$ for two adjacent subdomains $\Omega_{i_1}$ and $\Omega_{i_2}$

We consider the Poisson equation

$$
\begin{align}
  - \nabla^2 u = f \quad \text{in~$\Omega$}, \quad \quad 
  u = 0            \quad \text{on~$\partial \Omega$}, \quad \quad 
   [\![ u ]\!] = 0 \quad \text{on~$\mathcal{I}$}, \quad \quad
   [\![\partial_n u ]\!] = 0 \quad \text{on~$\mathcal{I}$}.
\end{align}
$$

The variational formulation reads

$$
\begin{align}
  \text{find $u \in V$ such that} \quad a(u,v) = l(v) \quad \forall v \in V,
\end{align}
$$

where 

- $V \subset H^1(\Omega)$, 
- $a(u,v) := \int_{\Omega} \nabla u \cdot \nabla v ~ d\Omega + \int_{\mathcal{I}} \left( \kappa [\![ u ]\!] ~ [\![ v ]\!] - [\![ u ]\!] ~ \partial_n v - \partial_n u ~ [\![ v ]\!] \right) ~ d\mathcal{I} $, 
- $l(v) := \int_{\Omega} f v ~ d\Omega ~ d\Gamma$.

## Examples

- [Poisson on a square with two subdomains and Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_poisson_subdomains_nitsche_dir0.py)


