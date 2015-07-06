M Model Collection
==================

A collection of M models downloaded from published studies, along with
notebooks which parse and solve them. The notebooks can be viewed at
[nbviewer](http://nbviewer.ipython.org/github/opencobra/m_model_collection/tree/master/).

To parse these reconstructions into a computing M-models, [install
cobrapy](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md)
version 0.4.0b1 or later, along with the optional dependcies libsbml and scipy.
To load the reconstructions from the Microsoft Excel format, the
[pandas](http://pandas.pydata.org/) library is also necessary. Then run
the ```load_models.py``` script.

The resulting ```all_models.mat``` file can be used with both cobrapy and the
COBRA toolbox. This script also generates SBML files which use the SBML level 3
with the fbc package (version 2), which are placed in the "sbml3" directory.

As shown in the included notebooks, these models can all be used to compute
[rational solutions](exact_solving_models.ipynb)
as well as [loopless solutions](loopless_m_models.ipynb).
They also can be solved by the [MATLAB cobra toolbox](solve_in_matlab_cobra_toolbox.ipynb).

Some of the [models](failed/README.md), however, could not be parsed correctly into
solving models. Still, they give [consistent results between floating point and
exact solvers](failed_models.ipynb)
