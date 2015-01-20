
# coding: utf-8

## Load and Process SBML models

# This script will load the M models in the collection using [libSBML](http://sbml.org/Software/libSBML) through cobrapy, and convert them to the "mat" format used by the COBRA toolbox (which can also be read and written by cobrapy).

# In[1]:

from os import listdir
import warnings
import re

import sympy
import scipy
import scipy.io

import cobra


# In addition to the usual fields in the "mat" struct, we will also include S_num and S_denom, which are the numerator and denominator of the stoichiometric coefficients encoded as rational numbers.

# In[2]:

def convert_to_rational(value):
    return sympy.Rational("%.15g" % value)


def construct_S_num_denom(model):
    """convert model to two S matrices

    they encode the numerator and denominator of stoichiometric
    coefficients encoded as rational numbers

    """
    # intialize to 0
    dimensions = (len(model.metabolites), len(model.reactions))
    S_num = scipy.sparse.lil_matrix(dimensions)
    S_denom = scipy.sparse.lil_matrix(dimensions)
    # populate with stoichiometry
    for i, r in enumerate(model.reactions):
        for met, value in r._metabolites.iteritems():
            rational_value = convert_to_rational(value)
            num, denom = (rational_value.p, rational_value.q)
            S_num[model.metabolites.index(met), i] = num
            S_denom[model.metabolites.index(met), i] = denom
    return S_num, S_denom


# There is quite a bit of code below to attempt to automatically identify model objectives when none is set by searching for reactions with "biomass" in a reaction or metabolite id. For some models however, the objective had to be determined manually. Additionally, some of the models need their exchange reactions opened.

# In[3]:

curated_objectives = {"VvuMBEL943": "R806",
                      "iAI549": "BIO_CBDB1_DM_855",
                      "mus_musculus": "BIO028",
                      "iRsp1095": "RXN1391",
                      "iLC915": "r1133",
                      "PpaMBEL1254": "R01288",
                      "AbyMBEL891": "R761"}
open_boundaries = {"iRsp1095", "AORYZAE_COBRA", "iFF708"}


# In[4]:

models = []
biomass_re = re.compile("biomass", re.IGNORECASE)
for i in sorted(listdir(".")):
    if not i.endswith(".xml"):
        continue
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = cobra.io.read_legacy_sbml(i)
    m.id = m.description = i[:-4]
    # Attempt to detect a biomass function when the model defines none
    if len(m.reactions.query(lambda x: x > 0, "objective_coefficient")) == 0:
        possible_objectives = m.reactions.query(biomass_re)
        # In some cases, a biomass "metabolite" is produced, whose production
        # should be the objective function.
        possible_biomass_metabolites = m.metabolites.query(biomass_re)
        if m.id in curated_objectives:
            m.change_objective(curated_objectives[m.id])
        elif len(possible_objectives) > 0:
            print("autodetected objective reaction '%s' for model '%s'" %
                  (possible_objectives[0].id, m.id))
            m.change_objective(possible_objectives[0])
        elif len(possible_biomass_metabolites) == 1:
            # In the case of a biomass metabolite, add a sink reaction for
            # it and make that the objective.
            biomass_met = possible_biomass_metabolites[0]
            r = cobra.Reaction("added_biomass_sink")
            r.objective_coefficient = 1
            r.add_metabolites({biomass_met: -1})
            m.add_reaction(r)
            print("autodetected biomass metabolite '%s' for model '%s'"
                  % (biomass_met.id, m.id))
        else:
            print("no objective found for " + m.id)
            continue
    # Ensure the biomass objective flux is unconstrained
    for reaction in m.reactions.query(lambda x: x > 0, "objective_coefficient"):
        reaction.lower_bound = min(reaction.lower_bound, 0)
        reaction.upper_bound = max(reaction.upper_bound, 1000)
    if m.id in open_boundaries:
        for reaction in m.reactions:
            if len(reaction.metabolites) == 1:
                # Ensure we are not creating any new sinks
                if reaction.metabolites.values()[0] > 0:
                    reaction.upper_bound = max(reaction.upper_bound, 10)
                else:
                    reaction.lower_bound = min(reaction.lower_bound, -10)
    models.append(m)


# Some models have erroneous formulas which are removed.

# In[5]:

for model in models:
    for metabolite in model.metabolites:
        if len(metabolite.formula.elements) == 0:
            metabolite.formula = ''


# Save all the models into a single mat file.

# In[6]:

all_model_dict = {}
for model in models:
    model_dict = cobra.io.mat.create_mat_dict(model)
    model_dict["S_num"], model_dict["S_denom"] = construct_S_num_denom(model)
    all_model_dict[model.id] = model_dict
scipy.io.savemat("all_models.mat", all_model_dict, oned_as="column")

