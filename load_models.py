
# coding: utf-8

# # Load and Process models
# 
# This script will load the M models in the collection using cobrapy, and convert them to a normalized format. They will also be exported to the "mat" format used by the COBRA toolbox.

# In[1]:

import os
import warnings
import re
from itertools import chain

import sympy
import scipy
import scipy.io

import cobra

from read_excel import read_excel


# ## Read in Models

# In[2]:

def open_exchanges(model, amount=10):
    for reaction in model.reactions:
        if len(reaction.metabolites) == 1:
            # Ensure we are not creating any new sinks
            if reaction.metabolites.values()[0] > 0:
                reaction.upper_bound = max(reaction.upper_bound, amount)
            else:
                reaction.lower_bound = min(reaction.lower_bound, -amount)

def add_exchanges(model, extracellular_suffix="[e]", uptake_amount=10):
    for metabolite in model.metabolites:
        if str(metabolite).endswith(extracellular_suffix):
            if len(metabolite.reactions) == 0:
                print "no reactions for " + metabolite.id
                continue
            if min(len(i.metabolites) for i in metabolite.reactions) > 1:
                EX_reaction = cobra.Reaction("EX_" + metabolite.id)
                EX_reaction.add_metabolites({metabolite: 1})
                m.add_reaction(EX_reaction)
                EX_reaction.upper_bound = uptake_amount
                EX_reaction.lower_bound = -uptake_amount


# ### SBML models
# 
# These models will be read in using [libSBML](http://sbml.org/Software/libSBML) through cobrapy.
# 
# There is quite a bit of code below to attempt to automatically identify model objectives when none is set by searching for reactions with "biomass" in a reaction or metabolite id. For some models however, the objective had to be determined manually. Additionally, some of the models need their exchange reactions opened.

# In[3]:

curated_objectives = {"VvuMBEL943": "R806",
                      "iAI549": "BIO_CBDB1_DM_855",
                      "mus_musculus": "BIO028",
                      "iRsp1095": "RXN1391",
                      "iLC915": "r1133",
                      "PpaMBEL1254": "R01288",
                      "AbyMBEL891": "R761",
                      "iAbaylyiV4": "GROWTH_DASH_RXN"}
open_boundaries = {"iRsp1095", "AORYZAE_COBRA", "iFF708"}
legacy_SBML = {"T_Maritima", "iNJ661m", "iSR432", "iTH366"}


# In[4]:

models = cobra.DictList()
biomass_re = re.compile("biomass", re.IGNORECASE)
for i in sorted(os.listdir("sbml")):
    if not i.endswith(".xml"):
        continue
    model_id = i[:-4]
    filepath = os.path.join("sbml", i)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = cobra.io.read_legacy_sbml(filepath) if model_id in legacy_SBML             else cobra.io.read_sbml_model(filepath)
    m.id = m.description = model_id.replace(".", "_")
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
        open_exchanges(m)
    models.append(m)


# ### Some models are only available as Microsoft Excel files

# In[5]:

m = read_excel("xls/iJS747.xls",
               verbose=False, rxn_sheet_header=7)
m.change_objective("agg_GS13m_2")
models.append(m)


# In[6]:

m = read_excel("xls/iRM588.xls",
               verbose=False, rxn_sheet_header=5)
m.change_objective("agg_GS13m")
models.append(m)


# In[7]:

m = read_excel("xls/iSO783.xls", verbose=False, rxn_sheet_header=2)
m.change_objective("Biomass")
models.append(m)


# In[8]:

m = read_excel("xls/iCR744.xls", rxn_sheet_header=4, verbose=False)
m.change_objective("BIO_Rfer3")
models.append(m)


# In[9]:

m = read_excel("xls/iNV213.xls", rxn_str_key="Reaction Formula", verbose=False)
m.change_objective("R_biomass_target")
# remove boundary metabolites
for met in list(m.metabolites):
    if met.id.endswith("[b]"):
        met.remove_from_model()
models.append(m)


# In[10]:

m = read_excel("xls/iTL885.xls", verbose=False,
               rxn_id_key="Rxn name", rxn_gpr_key="Gene-reaction association", met_sheet_name="ignore")
m.change_objective("SS1240")
models.append(m)


# In[11]:

m = read_excel("xls/iWZ663.xls", verbose=False,
               rxn_id_key="auto", rxn_name_key="Reaction name", rxn_gpr_key="Local gene")
m.change_objective(m.reactions.query("biomass", "name"))
models.append(m)


# In[12]:

m = read_excel("xls/iOR363.xls", verbose=False)
m.change_objective("OF14e_Retli")
models.append(m)


# In[13]:

m = read_excel("xls/iMA945.xls", verbose=False)
m.reactions.query("biomass")[0].objective_coefficient = 1.
models.append(m)


# In[14]:

m = read_excel("xls/iPP668.xls", verbose=False)
for met in list(m.metabolites):
    if met.id.endswith("[e]"):  # treats extracellular like boundary
        met.remove_from_model()
m.reactions.query("BIOMASS")[0].objective_coefficient = 1.
models.append(m)


# In[15]:

m = read_excel("xls/iVM679.xls", verbose=False, met_sheet_name="ignore",
               rxn_id_key="Name", rxn_name_key="Description", rxn_str_key="Reaction")
m.change_objective("Biomass")
open_exchanges(m)
models.append(m)


# In[16]:

m = read_excel("xls/iTY425.xls", rxn_sheet_header=1,
               rxn_sheet_name="S8", rxn_id_key="Number", rxn_str_key="Reaction", verbose=False)
m.change_objective("R517")
add_exchanges(m, "xt")
# PROTEIN component of biomass doesn't produce PROTEIN
m.reactions.R511.add_metabolites({m.metabolites.PROTEIN: 1})
# Remove BIOMASS metabolite
m.metabolites.BIOMASS.remove_from_model()
models.append(m)


# In[17]:

m = read_excel("xls/iSS724.xls", rxn_str_key="Reactions",
               rxn_sheet_header=1, met_sheet_header=1, rxn_id_key="Name",
               verbose=False)
m.change_objective("BIOMASS")
add_exchanges(m, "xt")
# Remove BIOMASS metabolite
m.metabolites.BIOMASS.remove_from_model()
models.append(m)


# ## Fixes of various encoding bugs

# ### General

# GSMN_TB does not use the convention of extracellular metabolites with exchanges. Although the model still solves with this formulation, this is still normalized here. This process does not change the mathematical structure of the model.

# In[18]:

h_c = models.GSMN_TB.metabolites.H_c
for r in models.GSMN_TB.reactions:
    if len(r.metabolites) == 2 and h_c in r.metabolites:
        met = [i for i in r.metabolites if i is not h_c][0]
        EX_met = cobra.Metabolite(met.id[:-1] + "e")
        r.add_metabolites({EX_met: -r.metabolites[met]})
        if "EX_" + EX_met.id not in models.GSMN_TB.reactions:
            exchange = cobra.Reaction("EX_" + EX_met.id)
            exchange.add_metabolites({EX_met: -1})
            exchange.lower_bound = -1000000.0
            exchange.upper_bound = 1000000.0
            models.GSMN_TB.add_reaction(exchange)


# ### Reaction and Metabolites

# ### id's

# In[19]:

# reaction id's with spaces in them
models.iJS747.reactions.get_by_id("HDH [deleted 01/16/2007  12:02:30 PM]").id = "HDH_del"
models.iJS747.reactions.get_by_id("HIBD [deleted 03/21/2007  01:06:12 PM]").id = "HIBD_del"
models.iAC560.reactions.get_by_id("GLUDx [m]").id = "GLUDx[m]"
for r in models.iOR363.reactions:
    if " " in r.id:
        r.id = r.id.split()[0]

models.textbook.reactions.query("Biomass")[0].id = "Biomass_Ecoli_core"


# Use the convention underscore + compartment i.e. _c instead of [c] (c) etc.

# In[20]:

SQBKT_re = re.compile("\[([a-z])\]$")

def fix_brackets(id_str, compiled_re):
    result = compiled_re.findall(id_str)
    if len(result) > 0:
        return compiled_re.sub("_" + result[0], id_str)
    else:
        return id_str

for r in models.iRS1597.reactions:
    r.id = fix_brackets(r.id, re.compile("_LSQBKT_([a-z])_RSQBKT_$"))

for m_id in ["iJS747", "iRM588", "iSO783", "iCR744", "iNV213", "iWZ663", "iOR363", "iMA945", "iPP668", "iTL885", "iVM679"]:
    for met in models.get_by_id(m_id).metabolites:
        met.id = fix_brackets(met.id, SQBKT_re)

for met in models.S_coilicolor_fixed.metabolites:
    if met.id.endswith("_None_"):
        met.id = met.id[:-6]

# Exchange reactions should have the id of the metabolite after with the same convention
for m_id in ["iAF1260", "iJO1366", "iAF692", "iJN746", "iRC1080", "textbook", "iNV213",
             "iIT341", "iJN678", "iJR904", "iND750", "iNJ661", "iPS189_fixed", "iSB619"]:
    for r in models.get_by_id(m_id).reactions:
        if len(r.metabolites) != 1:
            continue
        if r.id.startswith("EX_"):
            r.id = "EX_" + list(r.metabolites.keys())[0].id
        if r.id.startswith("DM_"):
            r.id = "DM_" + list(r.metabolites.keys())[0].id


# Ensure all id's are escaped

# In[21]:

def escape_id(id):
    id = id.replace("(", "_LPAREN_").replace(")", "_RPAREN_").replace("[", "_LSQBKT_").replace("]", "_RSQBKT_")
    id = id.replace(",", "__COMMA_").replace("/", "_FSLASH_").replace("_DASH_", "__").replace("-", "__")
    return id

for model in models:
    for x in chain(model.reactions, model.metabolites):
        x.id = escape_id(x.id)


# In[22]:

for m in models:
    m.repair()


# ### Metabolite Formulas

# In[23]:

for model in models:
    for metabolite in model.metabolites:
        if metabolite.formula is not None and str(metabolite.formula).lower() == "none":
            metabolite.formula = None
        # some characters should not be in a formula
        if metabolite.formula is not None:
            if "(" in metabolite.formula or                     ")" in metabolite.formula or                     "." in metabolite.formula:
                metabolite.formula = None


# ### Metabolite Compartments

# In[24]:

compartments = {
    'c': 'Cytoplasm',
    'e': 'Extracellular',
    'p': 'Periplasm',
    'm': 'Mitochondria',
    'g': 'Golgi',
    'n': "Nucleus",
    'r': "Endoplasmic reticulum",
    'x': "Peroxisome",
    'v': "Vacuole",
    "h": "Chloroplast",
    "x": "Glyoxysome",
    "s": "Eyespot",
    "default": "No Compartment"}

for model in models:
    for metabolite in model.metabolites:
        if metabolite.compartment is None or metabolite.compartment == "[":
            if len(metabolite.id) > 2 and metabolite.id[-2] == "_" and metabolite.id[-1].isalpha():
                metabolite.compartment = metabolite.id[-1]
            else:
                metabolite.compartment = "default"
        if metabolite.compartment not in model.compartments:
            model.compartments[metabolite.compartment] = compartments.get(metabolite.compartment, metabolite.compartment)


# ### Metabolite and Reaction Names
# Names which start with numbers don't need to be escaped with underscores.

# In[25]:

for model in models:
    for x in chain(model.metabolites, model.reactions):
        if x.name is not None and x.name.startswith("_"):
            x.name = x.name.lstrip("_")
        if x.name is not None:
            x.name = x.name.strip()


# ### MISC fixes

# In[26]:

models.iMM1415.reactions.EX_lnlc_dup_e.remove_from_model()
models.iMM1415.reactions.EX_retpalm_e.remove_from_model(remove_orphans=True)

# these reaction names are reaction strings
for r in models.iCac802.reactions:
    r.name = ""

for model_id in ["iRM588", "iJS747", "iCR744", "iSO783"]:
    model = models.get_by_id(model_id)
    for metabolite in model.metabolites:
        metabolite.id = metabolite.id.replace("4:2", "4_2")
    model.metabolites._generate_index()


# ## Fix Genes and GPR's

# A lot of genes have characters which won't work in their names

# In[27]:

def rename_gene(model, old_id, new_id, regexp=None):
    gene = model.genes.get_by_id(old_id)
    gene.id = new_id
    model.genes._dict[new_id] = model.genes._dict.pop(old_id)
    for reaction in gene.reactions:
        if regexp:
            new_gpr = regexp.sub(new_id, reaction._gene_reaction_rule)
        else:
            new_gpr = reaction._gene_reaction_rule.replace(old_id, new_id)
        reaction._gene_reaction_rule = new_gpr

# Fix N/A
rename_gene(models.iCac802, "N/A", "N_A")
rename_gene(models.iRS1563, "N/A", "N_A")

# nonbreaking spaces
for r in models.iCB925.reactions:
    if "\xa0" in r.gene_reaction_rule:
        r.gene_reaction_rule = r.gene_reaction_rule.replace("\xc2", " ").replace("\xa0", " ")
for g in list(models.iCB925.genes):
    if len(g.reactions) == 0:
        models.iCB925.genes.remove(g)

# dashes in ID's
for model_id in ["iLC915", "iND750", "iMM904", "iFF708", "AORYZAE_COBRA"]:
    model = models.get_by_id(model_id)
    for gene in model.genes:
        gene.name = gene.id
        gene.id = gene.id.replace("-", "__")
    model.repair()
    for reaction in model.reactions:
        reaction._gene_reaction_rule = reaction.gene_reaction_rule.replace("-", "__")

# just a number as gene id - needs to start with a character
for model_id in ["iAI549", "iMM1415"]:
    model = models.get_by_id(model_id)
    old_gene_ids = model.genes.list_attr("id")
    old_gene_ids.sort(key=lambda x: len(x))
    for old_id in old_gene_ids:
        # number should not be contained within another number
        regexp = re.compile("(?<!\d)%s(?!\d)" % old_id)
        rename_gene(model, old_id, "gene_" + old_id, regexp=regexp)


# Some GPR's have multiple ands/ors in a row

# In[28]:

multiple_ors = re.compile("(\s*or\s+){2,}")
multiple_ands = re.compile("(\s*and\s+){2,}")

for model_id in ["iRS1563", "iRS1597", "iMM1415"]:
    model = models.get_by_id(model_id)
    for reaction in model.reactions:
        gpr = reaction.gene_reaction_rule
        gpr = multiple_ors.sub(" or ", gpr)
        gpr = multiple_ands.sub(" and ", gpr)
        if "[" in gpr:
            gpr = gpr.replace("[", "(").replace("]", ")")
        if gpr.endswith(" or"):
            gpr = gpr[:-3]
        reaction.gene_reaction_rule = gpr
    for gene in list(model.genes):
        if gene.id.startswith("[") or gene.id.endswith("]"):
            if len(gene.reactions) == 0:
                model.genes.remove(gene.id)


# Some models are missing spaces between the ands/ors in some of their GPR's

# In[29]:

for m_id in ["iJN678", "iTL885"]:
    for r in models.get_by_id(m_id).reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.replace("and", " and ").replace("or", " or ").replace("  ", " ")


# Some models are missing parenthesis

# In[30]:

def insert(string, index, substr):
    return string[:index] + substr + string[index:]

models.iRS1597.reactions.R01138_c.gene_reaction_rule =     insert(models.iRS1563.reactions.R02235_c.gene_reaction_rule, 32, ")")
models.iRS1563.reactions.R02235_c.gene_reaction_rule =     insert(models.iRS1563.reactions.R02235_c.gene_reaction_rule, 32, ")")


# In[31]:

models.iCac802.reactions.R0095.gene_reaction_rule =     models.iCac802.reactions.R0095.gene_reaction_rule.replace(" AND ", " and ")


# ## Export Models

# ### mat
# Save all the models into a single mat file. In addition to the usual fields in the "mat" struct, we will also include S_num and S_denom, which are the numerator and denominator of the stoichiometric coefficients encoded as rational numbers.

# In[32]:

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


# In[33]:

all_model_dict = {}
for model in models:
    model_dict = cobra.io.mat.create_mat_dict(model)
    model_dict["S_num"], model_dict["S_denom"] = construct_S_num_denom(model)
    all_model_dict[model.id] = model_dict
scipy.io.savemat("all_models.mat", all_model_dict, oned_as="column")


# ### SBML 3
# Export the models to the use the draft fbc version 2 extension to SBML level 3 version 1.
# 
# Draft FBC 2 support is still [under development](https://github.com/opencobra/cobrapy/pull/152) in cobrapy.
from cobra.io.sbml3 import write_sbml_model

for model in models:
    write_sbml_model(model, "sbml3/%s.xml" % model.id)