# Analytical Mapping 


Analytical Mappings are provided as symbolic expressions, which allow us to compute automatically their jacobian matrices and all related geometrical information.

## Available mappings

| SymPDE objects | Physical Dimension | Logical Dimension | Description |
| -------------- | -----------------  | ----------------  | ----------- |
| `IdentityMapping` | 1D, 2D, 3D | 1D, 2D, 3D | Identity Mapping object <br> $\begin{cases} x&=x_1, \\ y&=x_2, \\ z &= x_3 \end{cases} $ |
| `AffineMapping`   | 1D, 2D, 3D | 1D, 2D, 3D | Affine Mapping object <br>  $\begin{cases} x &= c_1 + a_{11} x_1 + a_{12} x_2 + a_{13} x_3, \\ y &= c_2 + a_{21} x_1 + a_{22} x_2 + a_{23} x_3, \\ z &= c_3 + a_{31} x_1 + a_{32} x_2 + a_{33} x_3 \end{cases}$|
| `PolarMapping` | 2D | 2D | Polar Mapping object (Annulus) <br> $ \begin{cases} x &= c_1 + (r_{min} (1-x_1)+r_{max} x_1) \cos(x_2), \\ y &= c_2 + (r_{min} (1-x_1)+r_{max} x_1) \sin(x_2) \end{cases} $ |
| `TargetMapping` | 2D | 2D | Target Mapping object <br>  $\begin{cases} x &= c_1 + (1-k) x_1 \cos(x_2) - D x_1^2, \\ y &= c_2 + (1+k) x_1 \sin(x_2) \end{cases}$|
| `CzarnyMapping` | 2D | 2D | Czarny Mapping object <br>  $\begin{cases} x &= \frac{1}{\epsilon}(1 - \sqrt{ 1 + \epsilon (\epsilon + 2 x_1 \cos(x_2)) }),  \\  y &= c_2 + \frac{b}{\sqrt{1-\epsilon^2/4}} \frac{ x_1  \sin(x_2)}{2 - \sqrt{ 1 + \epsilon (\epsilon + 2 x_1 \cos(x_2)) }} \end{cases}$|
| `CollelaMapping2D` | 2D | 2D | Collela Mapping object <br>  $\begin{cases} x &= 2 (x_1 + \epsilon \sin(2 \pi k_1 x_1) \sin(2 \pi k_2 x_2)) - 1,  \\  y &= 2 (x_2 + \epsilon \sin(2 \pi k_1 x_1) \sin(2 \pi k_2 x_2)) - 1 \end{cases}$|
| `TorusMapping` | 3D | 3D | Parametrization of a torus (or a portion of it) <br>  $\begin{cases} x &= (R_0 + x_1 \cos(x_2)) \cos(x_3),  \\  y &= (R_0 + x_1 \cos(x_2)) \sin(x_3),  \\  z &= x_1   \sin(x_2) \end{cases}$|
| `TorusSurfaceMapping` | 3D | 2D  | surface obtained by "slicing" the torus above at r = a <br>  $\begin{cases} x &= (R_0 + a \cos(x_1)) \cos(x_2),   \\  y &= (R_0 + a \cos(x_1)) \sin(x_2),   \\  z &=       a \sin(x_1) \end{cases}$|
| `TwistedTargetSurfaceMapping` | 3D | 2D  | surface obtained by "twisting" the `TargetMapping` out of the (x, y) plane <br>  $\begin{cases} x &= c_1 + (1-k)  x_1  \cos(x_2) - D  x_1^2   \\  y &= c_2 + (1+k)  x_1  \sin(x_2)   \\  z &= c_3 + x_1^2 \sin(2 x_2) \end{cases}$|
| `TwistedTargetMapping` | 3D | 3D | volume obtained by "extruding" the `TwistedTargetSurfaceMapping` along z <br>  $\begin{cases} x &= c_1 + (1-k) x_1 \cos(x_2) - D x_1^2,   \\  y &= c_2 + (1+k) x_1 \sin(x_2),   \\  z &= c_3 + x_3  x_1^2 \sin(2 x_2) \end{cases}$|
| `SphericalMapping` | 3D | 3D | Parametrization of a sphere (or a portion of it) <br>  $\begin{cases} x &= x_1  \sin(x_2) \cos(x_3),   \\  y &= x_1  \sin(x_2) \sin(x_3),   \\  z &= x_1  \cos(x_2) \end{cases}$|

## Examples


```python
from sympde.topology      import Square
from sympde.topology      import PolarMapping

ldomain = Square('A',bounds1=(0., 1.), bounds2=(0, np.pi))
mapping = PolarMapping('M1',2, c1= 0., c2= 0., rmin = 0., rmax=1.)

domain = mapping(ldomain)

x,y = domain.coordinates
```


```python
mapping_1 = IdentityMapping('M1', 2)
mapping_2 = PolarMapping   ('M2', 2, c1 = 0., c2 = 0.5, rmin = 0., rmax=1.)
mapping_3 = AffineMapping  ('M3', 2, c1 = 0., c2 = np.pi, a11 = -1, a22 = -1, a21 = 0, a12 = 0)

A = Square('A',bounds1=(0.5, 1.), bounds2=(-1., 0.5))
B = Square('B',bounds1=(0.5, 1.), bounds2=(0, np.pi))
C = Square('C',bounds1=(0.5, 1.), bounds2=(np.pi-0.5, np.pi + 1))

D1     = mapping_1(A)
D2     = mapping_2(B)
D3     = mapping_3(C)

connectivity = [((0,1,1),(1,1,-1)), ((1,1,1),(2,1,-1))]
patches = [D1, D2, D3]
domain = Domain.join(patches, connectivity, 'domain')
```

```python
A = Square('A',bounds1=(0.2, 0.6), bounds2=(0, np.pi))
B = Square('B',bounds1=(0.2, 0.6), bounds2=(np.pi, 2*np.pi))
C = Square('C',bounds1=(0.6, 1.), bounds2=(0, np.pi))
D = Square('D',bounds1=(0.6, 1.), bounds2=(np.pi, 2*np.pi))

mapping_1 = PolarMapping('M1',2, c1= 0., c2= 0., rmin = 0., rmax=1.)
mapping_2 = PolarMapping('M2',2, c1= 0., c2= 0., rmin = 0., rmax=1.)
mapping_3 = PolarMapping('M3',2, c1= 0., c2= 0., rmin = 0., rmax=1.)
mapping_4 = PolarMapping('M4',2, c1= 0., c2= 0., rmin = 0., rmax=1.)

D1     = mapping_1(A)
D2     = mapping_2(B)
D3     = mapping_3(C)
D4     = mapping_4(D)

connectivity = [((0,1,1),(1,1,-1)), ((2,1,1),(3,1,-1)), ((0,0,1),(2,0,-1)),((1,0,1),(3,0,-1))]
patches = [D1, D2, D3, D4]
domain = Domain.join(patches, connectivity, 'domain')
```
