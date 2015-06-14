M Model Collection
==================

A collection of M models downloaded from published studies.

To parse these reconstructions into a computing M-models, [install
cobrapy](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md)
version 0.3.2 or later, along with the optional dependcies libsbml and scipy.
To load the reconstructions from the Microsoft Excel format, the
[pandas](http://pandas.pydata.org/) library is also necessary. Then run
the ```load_models.py``` script.

The resulting ```all_models.mat``` file can be used with both cobrapy and the
COBRA toolbox.

As shown in the included notebooks, these models can all be used to compute
[rational solutions](http://nbviewer.ipython.org/github/opencobra/m_model_collection/blob/master/exact_solving_models.ipynb)
as well as [loopless solutions](http://nbviewer.ipython.org/github/opencobra/m_model_collection/blob/master/loopless_m_models.ipynb).
They also can be solved by the [MATLAB cobra toolbox](http://nbviewer.ipython.org/github/opencobra/m_model_collection/blob/master/solve_in_matlab_cobra_toolbox.ipynb).
