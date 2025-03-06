# Boundary Conditions


SymPDE & Psydac allows you to use both strong and weak boundary conditions.
We start first by explaining how to identify a boundary in **NCube** domains, such as **Line**, **Square** and a **Cube**.

## Identifying a boundary in **NCube** domains

Since we're dealing with **NCube** domains, the best way to identify a boundary is to use the couple **(axis, extremity)**. The axis is usually defined as perpendicular to the boundary, while the **extremity** gives the lower or upper bound of the interval associated to the **axis**. 

Let's see it through the following examples;

For a **Line**, the following figure gives the value to access a boundary.

![png](images/boundary-conditions/line-boundaries.png)

For a **Square**, we have,

![png](images/boundary-conditions/square-boundaries.png)

For a **Cube**, it is similar and would not be helpful to visualize it here.

## Essential boundary conditions

Essential boundary conditions are usually treated in two ways, either in the strong or weak form. The latter can be achieved through the use of Nitsche's method. You can check the associated examples in the sequel.

For the essential boundary conditions, usually, one needs one of the following expressions;

$$
\begin{align}
u &= 0, \quad \text{(homogeneous case)}, \\
u &= f, \quad \text{(inhomogeneous case)}, \\
\mathbf{v} &= 0, \quad \text{(homogeneous case)}, \\
\mathbf{v} &= \mathbf{g}, \quad \text{(homogeneous case)}, 
\end{align}
$$

where $u$ and $f$ are scalar functions while $\mathbf{v}$ and $\mathbf{g}$ denote vector functions. 


Let `Gamma` denotes the boundary $\Gamma$.
The normal derivative is an available operator in SymPDE, it is provided as `Dn` operator. Alternatively, you can also define it manually,

```python
    nn = NormalVector('nn')
    dn = lambda a: dot(grad(a), nn)
```


we have the following equivalences,

| Mathematical Expression                                   | SymPDE Expression              |                Example             |
| --------------------------------------------------------- | ------------------------------ | ---------------------------------- |
| $u = 0$ on $\Gamma$                                       | `EssentialBC(u, 0, Gamma)`     |                                    |
| $u = f$ on $\Gamma$                                       | `EssentialBC(u, f, Gamma)`     |                                    |
| $\partial_n u = 0$ on $\Gamma$                            | `EssentialBC(dn(u), 0, Gamma)` |                                    |
| $\partial_n u = f$ on $\Gamma$                            | `EssentialBC(dn(u), f, Gamma)` |                                    |
| $\mathbf{v} = 0$ on $\Gamma$                              | `EssentialBC(v, 0, Gamma)`     |                                    |
| $\mathbf{v} = \mathbf{g}$ on $\Gamma$                     | `EssentialBC(v, g, Gamma)`     |                                    |
| $\mathbf{v} \cdot \mathbf{n} = 0$ on $\Gamma$             | ?                              |                                    |
| $\mathbf{v} \times \mathbf{n} = 0$ on $\Gamma$            | ?                              |                                    |
| $\mathbf{v} \cdot \mathbf{n} = f,\mathbf{g}$ on $\Gamma$  | ?                              |                                    |
| $\mathbf{v} \times \mathbf{n} = f,\mathbf{g}$ on $\Gamma$ | ?                              |                                    |

