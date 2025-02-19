# 1D Burgers equation


We consider the 1d Burgers equation

$$
\partial_t u + u \partial_x u = \nu \frac{\partial ^2u}{\partial x^2}
$$

$u_0(x) := u(x,t)$ denotes the initial condition.
We choose homogeneous neumann boundary conditions in this example, i.e.
$$\partial_n u = 0, \partial \Omega$$ with $\Omega = (0,1)$

## Time scheme
We shall use a $\theta$-scheme in this case and consider the following problem

$$
\begin{align}
\frac{u^{n+1}-u^n}{\Delta t} + 
\theta~ u^{n+1} \partial_x u^{n+1} + (1-\theta)~ u^n \partial_x u^n = \theta~\nu \frac{\partial ^2u^{n+1}}{\partial x^2} + (1-\theta)~\nu \frac{\partial ^2u^{n}}{\partial x^2}
\end{align}
$$

hence

$$
\begin{align}
u^{n+1} + \Delta t ~ \theta~ u^{n+1} \partial_x u^{n+1} - \Delta t ~ \theta~\nu \frac{\partial ^2u^{n+1}}{\partial x^2} =
u^{n} - \Delta t ~ (1-\theta)~ u^{n} \partial_x u^{n} + \Delta t ~ (1-\theta)~\nu \frac{\partial ^2u^{n}}{\partial x^2}
\end{align}
$$

from now on, we shall denote by $f^n$ the right hand side of the previous equation

$$f^n := u^{n} - \Delta t ~ (1-\theta)~ u^{n} \partial_x u^{n} + \Delta t ~ (1-\theta)~\nu \frac{\partial ^2u^{n}}{\partial x^2}$$

## Weak formulation

Let $v \in \mathcal{V}$ be a function test, we have by integrating by parts the highest order term:

$$
\begin{align}
\langle v, u^{n+1}  \rangle + \Delta t ~ \theta~ \langle v, u^{n+1} \partial_x u^{n+1}  \rangle + \Delta t ~ \theta~\nu \langle \frac{\partial v}{\partial x}, \frac{\partial u^{n+1}}{\partial x}  \rangle = \langle v, f^n \rangle
\end{align}
$$

The previous weak formulation is still nonlinear with respect to $u^{n+1}$. We shall then follow the same strategy as for the previous chapter on nonlinear Poisson problem.

The strategy is to define the left hand side as a **LinearForm** with respect to $v$, then linearize it around $u^{n+1}$. We therefor can use either Picard or Newton method to treat the nonlinearity.

We consider the following linear form

$$
\begin{align}
G(v;u,w) := \langle v, u  \rangle + \Delta t ~ \theta~ \langle v, w \partial_x u  \rangle + \Delta t ~ \theta~\nu \langle \frac{\partial v}{\partial x}, \frac{\partial u}{\partial x}  \rangle , \quad \forall u,v,w \in \mathcal{V}
\end{align}
$$

Our problem is then

$$
\begin{align}
\mbox{Find } u^{n+1} \in \mathcal{V}, \mbox{such that}\\
G(v;u^{n+1},u^{n+1}) = l(v), \quad \forall v \in \mathcal{V}
\end{align}
$$

where

$$
l(v) := \int_{\Omega} f^n v ~d\Omega, \quad \forall v \in \mathcal{V}
$$

## SymPDE code

```python
domain = Line()

V = ScalarFunctionSpace('V', domain)
u  = element_of(V, name='u')
v  = element_of(V, name='v')
w  = element_of(V, name='w')
un  = element_of(V, name='un') # time iteration
uk  = element_of(V, name='uk') # nonlinear solver iteration

x = domain.coordinates

nu = Constant('nu')
theta = Constant('theta')
dt = Constant('dt')
```

### Defining the Linear form $G$

```python
# Linear form g: V --> R
expr = v * u + dt*theta*v*w*dx(u) + dt*theta*nu*dx(v)*dx(u)
g = LinearForm(v, integral(domain, expr))
```

### Defining the Linear form $l$

```python
# Linear form l: V --> R
expr = v * un - dt*theta*v*un*dx(un) - dt*theta*nu*dx(v)*dx(un)
l = LinearForm(v, integral(domain, expr))
```


### Picard Method

$$
\begin{align}
\mbox{Find } u_{n+1} \in \mathcal{V}\_h, \mbox{such that}\\
G(v;u_{n+1},u_n) = l(v), \quad \forall v \in \mathcal{V}\_h
\end{align}
$$

#### Picard iteration

```python
# Variational problem
picard = find(u, forall=v, lhs=g(v, u=u,w=uk), rhs=l(v))
```

### Newton Method

Let's define 

$$
F(v;u) := G(v;u,u) -l(v), \quad \forall v \in \mathcal{V}
$$

Newton method writes

$$
\begin{align}
\mbox{Find } u_{n+1} \in \mathcal{V}\_h, \mbox{such that}\\
F^{\prime}(\delta u,v; u_n) = - F(v;u_n), \quad \forall v \in \mathcal{V} \\
u_{n+1} := u_{n} + \delta u, \quad \delta u \in \mathcal{V}
\end{align}
$$

#### Newton iteration

```python
F = LinearForm(v, g(v,w=u)-l(v))
du  = element_of(V, name='du')

Fprime = linearize(F, u, trials=du)

# Variational problem
newton = find(du, forall=v, lhs=Fprime(du, v,u=uk), rhs=-F(v,u=uk))
```

