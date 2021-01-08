# IGA-Python

This project provides a sequence of Jupyter notebooks introducing and using the **Isogeometric Analysis (IGA)** approach.

## Contents

### Part I - Computer Aided Design

[Introduction to CAD](https://nbviewer.jupyter.org/github/UM6P/Introduction-to-CAD/blob/main/README.md)

### Part II - Isogeometric Analysis

[Chapter 0 - Introduction](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter0/README.md)

[Chapter 1 - Linear problems](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter1/README.md)

[Chapter 2 - Non linear problems](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter2/README.md)

[Chapter 3 - Boundary conditions in IGA](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter3/README.md)

[Chapter 4 - Compatible Finite Elements](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter4/README.md)

[Chapter 5 - Navier-Stokes](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter5/README.md)

[Chapter 6 - Stabilized Finite Elements](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter6/README.md)

[Chapter 7 - Spectral Analysis using GLT](https://github.com/ratnania/IGA-Python/blob/main/lessons/Chapter7/README.md)

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

