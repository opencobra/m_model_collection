
# coding: utf-8

# # Load and Process models
# 
# This script will load the M models in the collection using cobrapy, and convert them to a normalized format. They will also be exported to the "mat" format used by the COBRA toolbox.
# 
# This requires [cobrapy](https://opencobra.github.io/cobrapy) version 0.4.0b1 or later.

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
# These models will be read in using [libSBML](http://sbml.org/Software/libSBML) through cobrapy. Some models will need their exchanges opened.

# In[3]:

legacy_SBML = {"T_Maritima", "iNJ661m", "iSR432", "iTH366"}
open_boundaries = {"iRsp1095", "iWV1314", "iFF708", "iZM363"}

models = cobra.DictList()
for i in sorted(os.listdir("sbml")):
    if not i.endswith(".xml"):
        continue
    model_id = i[:-4]
    filepath = os.path.join("sbml", i)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = cobra.io.read_legacy_sbml(filepath) if model_id in legacy_SBML             else cobra.io.read_sbml_model(filepath)
    m.id = m.description = model_id.replace(".", "_")

    if m.id in open_boundaries:
        open_exchanges(m)

    models.append(m)


# ### Models available in COBRA Toolbox "mat" format

# In[4]:

for i in sorted(os.listdir("mat")):
    if not i.endswith(".mat"):
        continue
    m = cobra.io.load_matlab_model(os.path.join("mat", i))
    m.id = i[:-4]
    
    if m.id in open_boundaries:
        open_exchanges(m)

    models.append(m)


# ### Some models are only available as Microsoft Excel files

# In[5]:

m = read_excel("xls/iJS747.xls",
               verbose=False, rxn_sheet_header=7)
models.append(m)


# In[6]:

m = read_excel("xls/iRM588.xls",
               verbose=False, rxn_sheet_header=5)
models.append(m)


# In[7]:

m = read_excel("xls/iSO783.xls", verbose=False, rxn_sheet_header=2)
models.append(m)


# In[8]:

m = read_excel("xls/iCR744.xls", rxn_sheet_header=4, verbose=False)
models.append(m)


# In[9]:

m = read_excel("xls/iNV213.xls", rxn_str_key="Reaction Formula", verbose=False)
# remove boundary metabolites
for met in list(m.metabolites):
    if met.id.endswith("[b]"):
        met.remove_from_model()
models.append(m)


# In[10]:

m = read_excel("xls/iTL885.xls", verbose=False,
               rxn_id_key="Rxn name", rxn_gpr_key="Gene-reaction association", met_sheet_name="ignore")
models.append(m)


# In[11]:

m = read_excel("xls/iWZ663.xls", verbose=False,
               rxn_id_key="auto", rxn_name_key="Reaction name", rxn_gpr_key="Local gene")
models.append(m)


# In[12]:

m = read_excel("xls/iOR363.xls", verbose=False)
models.append(m)


# In[13]:

m = read_excel("xls/iMA945.xls", verbose=False)
models.append(m)


# In[14]:

m = read_excel("xls/iPP668.xls", verbose=False)
add_exchanges(m)
models.append(m)


# In[15]:

m = read_excel("xls/iVM679.xls", verbose=False, met_sheet_name="ignore",
               rxn_id_key="Name", rxn_name_key="Description", rxn_str_key="Reaction")
open_exchanges(m)
models.append(m)


# In[16]:

m = read_excel("xls/iTY425.xls", rxn_sheet_header=1,
               rxn_sheet_name="S8", rxn_id_key="Number", rxn_str_key="Reaction", verbose=False)
add_exchanges(m, "xt")
# Protein production reaction does not prdocue "PROTEIN" metabolite
m.reactions.R511.add_metabolites({m.metabolites.PROTEIN: 1})
m.id = m.id + "_fixed"
models.append(m)


# In[17]:

m = read_excel("xls/iSS724.xls", rxn_str_key="Reactions",
               rxn_sheet_header=1, met_sheet_header=1, rxn_id_key="Name",
               verbose=False)
add_exchanges(m, "xt")
models.append(m)


# In[18]:

m = read_excel("xls/iCS400.xls", rxn_sheet_name="Complete Rxn List",
               rxn_sheet_header=2, rxn_str_key="Reaction",
               rxn_id_key="Name", verbose=False)
add_exchanges(m, "xt")
models.append(m)


# In[19]:

m = read_excel("xls/iLL672.xls",
               rxn_id_key="auto", met_sheet_name="Appendix 3 iLL672 metabolites",\
               rxn_str_key="REACTION", rxn_gpr_key="skip", verbose=False,
               rxn_sheet_name='Appendix 3 iLL672 reactions')
m.reactions[-1].objective_coefficient = 1
m.metabolites.BM.remove_from_model()
add_exchanges(m, "xt")
models.append(m)


# In[20]:

plus_re = re.compile("(?<=\S)\+")  # substitute H+ with H, etc.
m = read_excel("xls/iMH551.xls", rxn_sheet_name="GPR Annotation", rxn_sheet_header=4,
               rxn_id_key="auto", rxn_str_key="REACTION", rxn_gpr_key="skip",
               rxn_name_key="ENZYME", rxn_skip_rows=[625, 782, 787], verbose=False,
               rxn_sheet_converters={"REACTION": lambda x: plus_re.sub("", x)})
for met in m.metabolites:
    if met.id.endswith("(extracellular)"):
        met.id = met.id[:-15] + "_e"
m.repair()
add_exchanges(m, "_e")
models.append(m)


# In[21]:

m = read_excel("xls/iCS291.xls", rxn_sheet_name="Sheet1",
               rxn_str_key="Reaction",
               rxn_sheet_header=5, rxn_id_key="Name",
               verbose=False)
add_exchanges(m, "xt")
# BIOMASS is just all model metabolites in the Demands list
m.add_reaction(cobra.Reaction("BIOMASS"))
# taken from Table 1 in publication
biomass_mets = {}
for i in {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
          "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
          "THR", "TRP", "TYR", "VAL", "PTRC", "SPMD", "ATP", "GTP",
          "CTP", "UTP", "DATP", "DGTP", "DCTP", "DTTP", "PS", "PE",
          "PG", "PEPTIDO", "LPS", "OPP", "UDPP", "NAD", "NADP", "FAD",
          "COA", "ACP", "PTH", "THIAMIN", "MTHF", "MK", "DMK"
         }:
    biomass_mets[m.metabolites.get_by_id(i)] = -1
    dm = cobra.Reaction("DM_" + i)
    m.add_reaction(dm)
    dm.add_metabolites({m.metabolites.get_by_id(i): -1})
m.reactions.BIOMASS.add_metabolites(biomass_mets)
m.change_objective("BIOMASS")
add_exchanges(m, "xt")
models.append(m)


# In[22]:

m = read_excel("xls/iYO844.xls", rxn_sheet_name="Reaction and locus", verbose=False, rxn_gpr_key="Locus name",
           rxn_str_key=u'Equation (note [c] and [e] at the beginning refer to the compartment \n'
                       'the reaction takes place in, cytosolic and extracellular respectively)')
add_exchanges(m)

# create the biomass reaction from supplementary data table
# http://www.jbc.org/content/suppl/2007/06/29/M703759200.DC1/Biomass_composition.doc
r = cobra.Reaction("biomass")
r.objective_coefficient = 1.
m.add_reaction(r)
r.reaction = ("408.3 gly[c] + 266.9 ala-L[c] + 306.7 val-L[c] + 346.4 leu-L[c] + 269.9 ile-L[c] + "
              "216.2 ser-L[c] + 186.3 thr-L[c] + 175.9 phe-L[c] + 110.8 tyr-L[c] + 54.3 trp-L[c] + "
              "56.7 cys-L[c] + 113.3 met-L[c] + 323.1 lys-L[c] + 193.0 arg-L[c] + 81.7 his-L[c] + "
              "148.0 asp-L[c] +  260.4 glu-L[c] + 148.0 asp-L[c] + 260.3 gln-L[c] + 160.6 pro-L[c] + "
              "62.7 gtp[c] + 38.9 ctp[c] + 41.5 utp[c] + 23.0 datp[c] + 17.4 dgtp[c] + 17.4 dctp[c] + "
              "22.9 dttp[c] + 0.085750 m12dg_BS[c] + 0.110292 d12dg_BS[c] + 0.065833 t12dg_BS[c] + "
              "0.004642 cdlp_BS[c] + 0.175859 pgly_BS[c] + 0.022057 lysylpgly_BS[c] + 0.559509 psetha_BS[c] + "
              "0.006837 lipo1-24_BS[c] + 0.006123 lipo2-24_BS[c] + 0.018162 lipo3-24_BS[c] + "
              "0.014676 lipo4-24_BS[c] + 101.82 peptido_BS[c] + 3.62 gtca1-45_BS[c] + 2.35 gtca2-45_BS[c] + "
              "1.82 gtca3-45_BS[c] + 3.11 tcam_BS[c] + 706.3 k[c] + 101.7 mg2[c] + 3.4 fe3[c] + 3.2 ca2[c] + "
              "0.9 ppi[c] + 0.3 mql7[c] + 0.4 10fthf[c] + 16.2 nad[c] + 4.7 amp[c]  + 2.6 adp[c] + 1.0 cmp[c]  + "
              "0.9 nadp[c] + 0.5 ctp[c]  + 0.5 gmp[c] + 0.4 gtp[c] + 0.3 cdp[c] + 0.2 nadph[c] + 0.2 gdp[c] + "
              "105053.5 atp[c] + 105000 h2o[c] --> 104985.6 pi[c] + 104997.4 adp[c] + 105000 h[c]")
# units are in mg for this reaction, so scale to grams
r *= 0.001
models.append(m)


# In[23]:

models.sort()


# ## Determine Objective Reactions
# 
# Some of these models do not specify an objective (or biomass) reaction. These will be automatically detected if possible, or set from a manually curated list.

# In[24]:

# regular expression to detect "biomass"
biomass_re = re.compile("biomass", re.IGNORECASE)

# manually identified objective reactions
curated_objectives = {"VvuMBEL943": "R806",
                      "iAI549": "BIO_CBDB1_DM_855",
                      "mus_musculus": "BIO028",
                      "iRsp1095": "RXN1391",
                      "iLC915": "r1133",
                      "PpaMBEL1254": "R01288",
                      "AbyMBEL891": "R761",
                      "iAbaylyiV4": "GROWTH_DASH_RXN",
                      "iOG654": "RM00001",
                      "iOR363": "OF14e_Retli",
                      "iRM588": "agg_GS13m",
                      "iJS747": "agg_GS13m_2",
                      "iTL885": "SS1240",
                      "iMH551": "R0227"}


# In[25]:

for m in models:
    if len(m.reactions.query(lambda x: x > 0, "objective_coefficient")):
        continue
    if m.id in curated_objectives:
        m.change_objective(curated_objectives[m.id])
        continue
    
    # look for reactions with "biomass" in the id or name
    possible_objectives = m.reactions.query(biomass_re)
    if len(possible_objectives) == 0:
        possible_objectives = m.reactions.query(biomass_re, "name")
    
    # In some cases, a biomass "metabolite" is produced, whose production
    # should be the objective function.
    possible_biomass_metabolites = m.metabolites.query(biomass_re)
    if len(possible_biomass_metabolites) == 0:
        possible_biomass_metabolites = m.metabolites.query(biomass_re, "name")

    if len(possible_biomass_metabolites) > 0:
        biomass_met = possible_biomass_metabolites[0]
        r = cobra.Reaction("added_biomass_sink")
        r.objective_coefficient = 1
        r.add_metabolites({biomass_met: -1})
        m.add_reaction(r)
        print ("autodetected biomass metabolite '%s' for model '%s'" %
              (biomass_met.id, m.id))

    elif len(possible_objectives) > 0:
        print("autodetected objective reaction '%s' for model '%s'" %
              (possible_objectives[0].id, m.id))
        m.change_objective(possible_objectives[0])

    else:
        print("no objective found for " + m.id)

# Ensure the biomass objective flux is unconstrained
for m in models:
    for reaction in m.reactions.query(lambda x: x > 0, "objective_coefficient"):
        reaction.lower_bound = min(reaction.lower_bound, 0)
        reaction.upper_bound = max(reaction.upper_bound, 1000)


# ## Fixes of various encoding bugs

# ### General

# GSMN_TB does not use the convention of extracellular metabolites with exchanges. Although the model still solves with this formulation, this is still normalized here. This process does not change the mathematical structure of the model.

# In[26]:

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

# In[27]:

# reaction id's with spaces in them
models.iJS747.reactions.get_by_id("HDH [deleted 01/16/2007  12:02:30 PM]").id = "HDH_del"
models.iJS747.reactions.get_by_id("HIBD [deleted 03/21/2007  01:06:12 PM]").id = "HIBD_del"
models.iAC560.reactions.get_by_id("GLUDx [m]").id = "GLUDx[m]"
for r in models.iOR363.reactions:
    if " " in r.id:
        r.id = r.id.split()[0]

models.textbook.reactions.query("Biomass")[0].id = "Biomass_Ecoli_core"


# Use the convention underscore + compartment i.e. _c instead of [c] (c) etc.

# In[28]:

SQBKT_re = re.compile("\[([a-z])\]$")

def fix_brackets(id_str, compiled_re):
    result = compiled_re.findall(id_str)
    if len(result) > 0:
        return compiled_re.sub("_" + result[0], id_str)
    else:
        return id_str

for r in models.iRS1597.reactions:
    r.id = fix_brackets(r.id, re.compile("_LSQBKT_([a-z])_RSQBKT_$"))

for m_id in ["iJS747", "iRM588", "iSO783", "iCR744", "iNV213", "iWZ663", "iOR363", "iMA945", "iPP668",
             "iTL885", "iVM679", "iYO844", "iZM363"]:
    for met in models.get_by_id(m_id).metabolites:
        met.id = fix_brackets(met.id, SQBKT_re)

for met in models.S_coilicolor_fixed.metabolites:
    if met.id.endswith("_None_"):
        met.id = met.id[:-6]

# Some models only have intra and extracellular metabolites, but don't use _c and _e.
for m_id in ["iCS291", "iCS400", "iTY425_fixed", "iSS724"]:
    for metabolite in models.get_by_id(m_id).metabolites:
        if metabolite.id.endswith("xt"):
            metabolite.id = metabolite.id[:-2] + "_e"
        elif len(metabolite.id) < 2 or metabolite.id[-2] != "_":
            metabolite.id = metabolite.id + "_c"

# Exchange reactions should have the id of the metabolite after with the same convention
for m_id in ["iAF1260", "iJO1366", "iAF692", "iJN746", "iRC1080", "textbook", "iNV213",
             "iIT341", "iJN678", "iJR904", "iND750", "iNJ661", "iPS189_fixed", "iSB619",
             "iZM363", "iMH551"]:
    for r in models.get_by_id(m_id).reactions:
        if len(r.metabolites) != 1:
            continue
        if r.id.startswith("EX_"):
            r.id = "EX_" + list(r.metabolites.keys())[0].id
        if r.id.startswith("DM_"):
            r.id = "DM_" + list(r.metabolites.keys())[0].id

for m in models:
    m.repair()


# ### Metabolite Formulas

# In[29]:

for model in models:
    for metabolite in model.metabolites:
        if metabolite.formula is None:
            metabolite.formula = ""
            continue
        if str(metabolite.formula).lower() == "none":
            metabolite.formula = ""
            continue
        # some characters should not be in a formula
        if "(" in metabolite.formula or                 ")" in metabolite.formula or                 "." in metabolite.formula:
            metabolite.formula = ""


# ### Metabolite Compartments

# In[30]:

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
        if metabolite.compartment is None or len(metabolite.compartment.strip()) == 0 or metabolite.compartment == "[":
            if len(metabolite.id) > 2 and metabolite.id[-2] == "_" and metabolite.id[-1].isalpha():
                metabolite.compartment = metabolite.id[-1]
            else:
                metabolite.compartment = "default"
        if metabolite.compartment not in model.compartments:
            model.compartments[metabolite.compartment] = compartments.get(metabolite.compartment, metabolite.compartment)


# ### Metabolite and Reaction Names
# Names which start with numbers don't need to be escaped with underscores.

# In[31]:

for model in models:
    for x in chain(model.metabolites, model.reactions):
        if x.name is not None and x.name.startswith("_"):
            x.name = x.name.lstrip("_")
        if x.name is not None:
            x.name = x.name.strip()
        if x.name is None:
            x.name = x.id


# ### MISC fixes

# In[32]:

models.iMM1415.reactions.EX_lnlc_dup_e.remove_from_model()
models.iMM1415.reactions.EX_retpalm_e.remove_from_model(remove_orphans=True)

# these reaction names are reaction strings
for r in models.iCac802.reactions:
    r.name = ""


# ## Fix Genes and GPR's

# A lot of genes have characters which won't work in their names

# In[33]:

# nonbreaking spaces
models.iCB925.reactions.FDXNRy.gene_reaction_rule = '( Cbei_0661 or Cbei_2182 )'
for r in models.iCB925.reactions:
    if "\xa0" in r.gene_reaction_rule:
        r.gene_reaction_rule = r.gene_reaction_rule.replace("\xc2", " ").replace("\xa0", " ")
for g in list(models.iCB925.genes):
    if len(g.reactions) == 0:
        models.iCB925.genes.remove(g)


# Some GPR's are not valid boolean expressions.

# In[34]:

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
        if gpr.count("(") != gpr.count(")"):
            gpr = ""  # mismatched parenthesis somewhere
        reaction.gene_reaction_rule = gpr
    for gene in list(model.genes):
        if gene.id.startswith("[") or gene.id.endswith("]"):
            if len(gene.reactions) == 0:
                model.genes.remove(gene.id)
                
# Some models are missing spaces between the ands/ors in some of their GPR's
for m_id in ["iJN678", "iTL885"]:
    for r in models.get_by_id(m_id).reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.replace("and", " and ").replace("or", " or ")

models.iCac802.reactions.R0095.gene_reaction_rule =     models.iCac802.reactions.R0095.gene_reaction_rule.replace(" AND ", " and ")


# make sbml3 output deterministic by sorting genes
for m in models:
    m.genes.sort()


# ## Ensure all ID's are SBML compliant

# In[35]:

for m in models:
    cobra.manipulation.escape_ID(m)


# ## Export Models

# ### SBML 3
# Export the models to the use the fbc version 2 (draft RC6) extension to SBML level 3 version 1.

# In[36]:

for model in models:
    cobra.io.write_sbml_model(model, "sbml3/%s.xml" % model.id)


# ### mat
# Save all the models into a single mat file. In addition to the usual fields in the "mat" struct, we will also include S_num and S_denom, which are the numerator and denominator of the stoichiometric coefficients encoded as rational numbers.

# In[37]:

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


# In[38]:

all_model_dict = {}
for model in models:
    model_dict = cobra.io.mat.create_mat_dict(model)
    model_dict["S_num"], model_dict["S_denom"] = construct_S_num_denom(model)
    all_model_dict[model.id] = model_dict
scipy.io.savemat("all_models.mat", all_model_dict, oned_as="column")

