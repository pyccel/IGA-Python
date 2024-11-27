# Installation
*Author: Ahmed Ratnani*

It is recommanded to use a virtual environement.


```shell
# create a virtual environement
python3 -m venv iga-python-env

# activate this environement using 
source iga-python-env/bin/activate
```

Then install dependencies using

```shell
pip3 install 'psydac[extra] @ git+https://github.com/pyccel/psydac.git@devel#egg=psydac'
pip3 install notebook
```

Adding the virtual environement to Jupyter notebook, can be done using

```shell
python3 -m ipykernel install --user --name=.iga-python
```
