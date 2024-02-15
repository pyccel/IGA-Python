# Installation

It is recommanded to use a virtual environement.


```shell
# create a virtual environement
python3 -m venv .iga-python

# activate this environement using 
source .iga-python/bin/activate
```

Then install dependencies using

```shell
pip3 install wheel
# numba does work right now with numpy > 1.21
pip3 install "numpy<=1.21"
pip3 install pyccel 
pip3 install ipykernel
pip3 install git+https://github.com/pyccel/sympde.git
pip3 install git+https://github.com/pyccel/gelato.git
pip3 install git+https://github.com/pyccel/psydac.git
pip3 install git+https://github.com/girving/igakit.git
```

Adding the virtual environement to Jupyter notebook, can be done using

```shell
python3 -m ipykernel install --user --name=.iga-python
```
