M Model Collection
==================

To view the publications and organisms associated with each model, please see
the [key](model_key.md).

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

Unofficial Models
-----------------
Please note that these models are not the official versions. Many of the models
in here have been updated at the original source. This is merely a snapshot of
models which were used in the analysis for
[this study](http://dx.doi.org/10.15252/msb.20156157).

LICENSE
-------
All code in this repository is released into the public domain under the
unlicense. Please see http://unlicense.org for details. The work may also
be used under the terms of [CC0](http://creativecommons.org/publicdomain/zero/1.0/).

The models themselves, however, are provided only for reproducing the results
in the included scripts. All other uses must be consistent with the initial
license the model was released under.
