m_model_collection
==================

A collection of M models downloaded from published studies

To parse these SBML-encoded reconstructions into a computing M-models, please
[install cobrapy](https://github.com/opencobra/cobrapy/blob/master/INSTALL.md)
version 0.3.2 or later, along with the optional dependcies libsbml and scipy.
Then run the ```convert_models_to_mat.py``` script.

The resulting ```all_models.mat``` file can be used with both cobrapy and the
COBRA toolbox.

As shown in the included notebooks, these models can all be used to compute
[rational solutions](http://nbviewer.ipython.org/github/opencobra/m_model_collection/blob/master/exact_solving_models.ipynb)
as well as [loopless solutions](http://nbviewer.ipython.org/github/opencobra/m_model_collection/blob/master/loopless_m_models.ipynb).
