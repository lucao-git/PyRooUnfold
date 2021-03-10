# PyRooUnfold
This is a Python wrapper of RooUnfold, plus a toolkit for bias, toy test, and converter for dataframe etc. which are commonly used in Python. 



## Install

Please make sure RooUnfold has been installed before using PyRooUnfold.

If not, please first follow the steps to install RooUnfold https://gitlab.cern.ch/RooUnfold/RooUnfold/-/blob/master/README.md

After installing RooUnfold, please set the environment variable for your RooUnfold libary as

```
export ROOUNFOLD_PATH="/path/for/your/libRooUnfold.so"
```

e.g. in my case

```
export ROOUNFOLD_PATH="/Users/caolu/Workspace/RooUnfold/libRooUnfold.so"
```

To install PyRooUnfold, you can do

```
python3 setup.py install
```

If you do not have your own python installation you can do

```
python3 setup.py install --user
```


### Tipps on RooUnfold

Details of RooUnfold can be found at
https://gitlab.cern.ch/RooUnfold/RooUnfold and an older version at http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html


This project has been set up using PyScaffold 3.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
