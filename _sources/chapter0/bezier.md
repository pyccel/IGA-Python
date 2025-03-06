# Bernstein polynomials


Without loss of generality, we restrict to the case of the unit interval, namely $a=0$ and $b=1$.
In figure (Fig. \ref{fig:bernstein-polynomials}), we plot the first sixth Bernstein polynomials.


````{prf:definition}

Fro $n \in \mathbb{N}$, the $n$-th degree Bernstein polynomials $B_{j,n}\,:\, [0,1] \longrightarrow \mathbb{R}$ are defined by

$$
B_{j,n}(x)=\frac{n!}{(n-j)!j!} x^j (1-x)^{n-j}, \quad x \in [0,1], \quad j=0, \ldots,n.
$$
````

Bernstein polynomials exhibit interesting properties highlighted in the following proposition,

```python
# needed imports
import numpy as np
import matplotlib.pyplot as plt
```

We first consider the evaluation of the Bernstein polynomials. 

## Evaluation of Bernstein polynomials

The following function evaluates all Bernstein polynomials of degree $n$ at $x$

```python
def all_bernstein(n, x):
    b = np.zeros(n+1)
    b[0] = 1.
    x1 = 1.-x
    for j in range(1, n+1):
        saved = 0.
        for i in range(0, j):
            tmp = b[i]
            b[i] = saved + x1*tmp
            saved = x*tmp
        b[j] = saved
    return b
```

## Bernstein polynomials properties

````{prf:proposition}

We have the following properties of Bernstein Polynomials,

- Positivity: $B_{j,n}(x) \ge 0$, for all $x \in \left[ 0, 1 \right]$ and $0 \leq j \leq n$.
- Partition of unity: $\sum_{j=0}^n B_{j,n}(x) = 1$, for all $x \in \left[ 0, 1 \right]$.
- $B_{0,n}(0) = B_{n,n}(1) = 1$. 
- $B_{j,n}$ has exactly one maximum on the interval $\left[ 0, 1 \right]$, at $\frac{j}{n}$. 
- Symmetry: $B_{j,n}$ is symmetric with respect to $x = \frac{1}{2}$ for $j=0,\ldots,n$.
- Bernstein polynomials can be defined recursively using the formula 

$$
\begin{align}
  B_{j,n}(x) = (1-x) B_{j,n-1}(x) + x B_{k-1,n-1}(x),
\end{align}
$$

  where we assume $B_{j,n}(x) = 0$ if $j < 0$ or $j > n$.

- Bernstein derivatives can be computed using the formulae  

$$
\begin{align}
  {B_{j,n}}^\prime(x) = n \left(B_{j-1,n-1}(x) - B_{j,n-1}(x) \right).
\end{align}
$$
using the same assumption as before.
````
