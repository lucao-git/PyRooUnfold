# PyRooUnfold
This is a Python wrapper of RooUnfold, plus a toolkit for bias, toy test, and converter for dataframe etc. which are commonly used in Python. 



## Install

- Please make sure RooUnfold has been installed before using PyRooUnfold.

If not, please first follow the steps to install RooUnfold https://gitlab.cern.ch/RooUnfold/RooUnfold/-/blob/master/README.md

or, use a pre-installed RooUnfold (e.g. from CVMFS).  


- After installing RooUnfold, please set the environment variable for your RooUnfold libary as

```
export ROOUNFOLD_PATH="/path/for/your/libRooUnfold.so"
```

e.g. in my case

```
export ROOUNFOLD_PATH="/Users/caolu/Workspace/RooUnfold/libRooUnfold.so"
```

or using CVMFS version provided as externals of Belle-II software 

```
export ROOUNFOLD_PATH="/cvmfs/belle.cern.ch/el9/externals/v02-01-00/Linux_x86_64/common/lib/libRooUnfold.so"
```

- To install PyRooUnfold, you can do

```
pip install .
```

If stall via setup.py, you need to include the package into your python library

```
export PYTHONPATH=$PYTHONPATH:<path-of-PyrooUnfold>
```

### Versions
- v1.0.0: compatible with RooUnfold 2.0.1
- v2.0.0: compatible with RooUnfold 3.0.0 

#### Tipps on RooUnfold

Details of RooUnfold can be found at
https://gitlab.cern.ch/RooUnfold/RooUnfold and an older version at http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html






This project has been set up using PyScaffold 3.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
