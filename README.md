# IGA-Python

This project provides a sequence of Jupyter notebooks introducing and using the **Isogeometric Analysis (IGA)** approach.

## Contents

### Part I - Introduction to Computer Aided Design

[Chapter 1 - Bézier curves](https://nbviewer.jupyter.org/github/UM6P/Introduction-to-CAD/blob/main/notebooks/Bezier_curves.ipynb)

[Chapter 2 - B-Splines](https://nbviewer.jupyter.org/github/UM6P/Introduction-to-CAD/blob/main/notebooks/B-Splines.ipynb)

[Chapter 3 - B-Splines curves](https://nbviewer.jupyter.org/github/UM6P/Introduction-to-CAD/blob/main/notebooks/B-Splines_curves.ipynb)

[Chapter 4 - B-Splines Fundamental geometric operations](https://nbviewer.jupyter.org/github/UM6P/Introduction-to-CAD/blob/main/notebooks/B-Splines_Fundamental_geometric_operations.ipynb)

### Part II - Isogeometric Analysis

[Chapter 1 - Introduction](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter1/README.md)
          
[Chapter 2 - Linear problems](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter2/README.md)
          
[Chapter 3 - Non linear problems](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter3/README.md)
          
[Chapter 4 - Boundary conditions in IGA](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter4/README.md)
          
[Chapter 5 - Compatible Finite Elements](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter5/README.md)
          
[Chapter 6 - Navier-Stokes](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter6/README.md)
          
[Chapter 7 - Stabilized Finite Elements](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter7/README.md)

[Chapter 8 - Spectral Analysis using GLT](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter8/README.md)

## How to use this material

### Using a virtual environement (recommanded)

First you need to create a virtual environement:

```shell
python3 -m venv iga-python
```

you activate this environement using

```shell
source iga-python/bin/activate
```

Then install dependencies using

```shell
pip3 install -r requirements.txt
```

Adding the virtual environement to Jupyter notebook, can be done using

```shell
python3 -m ipykernel install --user --name=iga-python
```

### without virtual environement

Install dependencies using

```shell
pip3 install -r requirements.txt
```

