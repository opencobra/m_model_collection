from warnings import warn

from six import string_types
from math import isnan
import pandas
import os

from cobra import Model, Metabolite, Reaction

# Potential string keys in the excel files
RXN_SHEET_NAMES = {"reaction list", "reactions"}
MET_SHEET_NAMES = {"metabolite list", "metabolites"}
MET_ID_KEYS = {"abbreviation", "abbreviations", "metabolite abbreviation"}
MET_NAME_KEYS = {"metabolite name", "officialname", "metabolite_name"}
MET_FORMULA_KEYS = {"chemical formula", "formula"}
RXN_ID_KEYS = {"abbreviation", "reaction #", "abbrev",
               "reaction id", "reaction abbreviation"}
RXN_NAME_KEYS = {"name", "Rxn description"}
RXN_STR_KEYS = {"equation", "reaction formula"}
RXN_GPR_KEYS = {"gene", "geneassociation"}
RXN_LB_KEYS = {"lb", "lower bound", "lower bounds"}
RXN_UB_KEYS = {"ub", "upper bound", "upper bounds"}


def escape_str(potential_str):
    try:
        return str(potential_str)
    except UnicodeEncodeError:
        return potential_str.encode("ascii", "replace")


def guess_name(potential_names, allowed_names, fail=True):
    """see if any of the potential names are allowed for a category"""
    for potential_name in potential_names:
        if potential_name.lower() in allowed_names:
            return potential_name
    if fail:
        raise ValueError("could not find any of %s in %s" %
                         (str(set(allowed_names)),
                          str(set(potential_names))))
    else:
        return None


def extract(row, keyname, type=str):
    """extract a value which may be missing, for a potentially missing key"""
    if keyname is None:
        return type()
    value = row[keyname]
    return type() if isinstance(value, float) and isnan(value) else type(value)


def read_excel(
        filename,
        verbose=True,
        rxn_sheet_name=None,
        rxn_sheet_header=0,

        rxn_id_key=None,
        rxn_name_key=None,
        rxn_str_key=None,
        rxn_gpr_key=None,
        rxn_lb_key=None,
        rxn_ub_key=None,
        rxn_fwd_arrow=None,
        rxn_rev_arrow=None,
        rxn_reversible_arrow=None,

        met_sheet_name=None,
        met_id_key=None,
        met_name_key=None,
        met_formula_key=None,
        ):

    # autodetect sheet names
    pio = pandas.io.excel.ExcelFile(filename)
    sheet_names = pio.sheet_names
    if rxn_sheet_name is None:
        # only one sheet means it must be that one
        if len(sheet_names) == 1:
            rxn_sheet_name = sheet_names[0]
        else:
            rxn_sheet_name = guess_name(pio.sheet_names, RXN_SHEET_NAMES)
    if met_sheet_name is None:
        met_sheet_name = guess_name(
            pio.sheet_names, MET_SHEET_NAMES, fail=False)

    # Begin model creation
    try:
        model_id = os.path.splitext(os.path.split(filename)[1])[0]
    except:
        model_id = "imported_model"
    m = Model(model_id)

    # Metabolites - if the sheet is found
    if met_sheet_name is not None:
        met_frame = pandas.read_excel(filename, met_sheet_name)
        # strip spaces from header
        met_frame.columns = [i.strip() for i in met_frame.keys()]
        if met_id_key is None:
            met_id_key = guess_name(met_frame.keys(), MET_ID_KEYS)
        if met_name_key is None:
            met_name_key = guess_name(met_frame.keys(), MET_NAME_KEYS,
                                      fail=False)
        if met_formula_key is None:
            met_formula_key = guess_name(met_frame.keys(), MET_FORMULA_KEYS,
                                         fail=False)
        for i in met_frame.index:
            met_row = met_frame.ix[i]
            met_attributes = {}
            if met_formula_key is not None:
                formula = extract(met_row, met_formula_key)
                if formula is not None and formula.lower() != "None":
                    met_attributes["formula"] = formula
            met = Metabolite(met_row[met_id_key],
                             name=extract(met_row, met_id_key),
                             **met_attributes)
            try:
                m.add_metabolites(met)
            except ValueError:
                if verbose:
                    print("duplicate metabolite '%s' not added" % met.id)
    elif verbose:
        print("metabolite sheet not found")

    # Reactions
    rxn_frame = pandas.read_excel(filename, rxn_sheet_name,
                                  header=rxn_sheet_header)
    # strip spaces from header
    rxn_frame.columns = [i.strip() for i in rxn_frame.keys()]
    if rxn_id_key is None:
        rxn_id_key = guess_name(rxn_frame.keys(), RXN_ID_KEYS)
    if rxn_str_key is None:
        rxn_str_key = guess_name(rxn_frame.keys(), RXN_STR_KEYS)
    if rxn_name_key is None:
        rxn_name_key = guess_name(rxn_frame.keys(), RXN_NAME_KEYS, fail=False)
        if verbose and rxn_name_key is None:
            print("reaction name column not identified")
    if rxn_gpr_key is None:
        rxn_gpr_key = guess_name(rxn_frame.keys(), RXN_GPR_KEYS, fail=False)
        if verbose and rxn_gpr_key is None:
            print("gene reaction rule column not identified")
    if rxn_lb_key is None:
        rxn_lb_key = guess_name(rxn_frame.keys(), RXN_LB_KEYS, fail=False)
        if verbose and rxn_lb_key is None:
            print("reaction lower bound column not identified")
    if rxn_ub_key is None:
        rxn_ub_key = guess_name(rxn_frame.keys(), RXN_UB_KEYS, fail=False)
        if verbose and rxn_ub_key is None:
            print("reaction upper bound column not identified")

    for i in range(len(rxn_frame)):
        row = rxn_frame.ix[i]
        rxn_id = row[rxn_id_key]
        rxn_str = row[rxn_str_key]
        if not isinstance(rxn_id, string_types) or \
                not isinstance(rxn_str, string_types):
            continue
        rxn = Reaction()
        rxn.id = escape_str(rxn_id)
        rxn.name = extract(row, rxn_name_key)
        if rxn.id in m.reactions:
            if verbose:
                print("duplicate reaction '%s' found" % rxn.id)
            # add underscores to the end of duplicate reactions
            while rxn.id in m.reactions:
                rxn.id += "_"
        m.add_reaction(rxn)
        rxn.build_reaction_from_string(
            rxn_str, verbose=verbose, fwd_arrow=rxn_fwd_arrow,
            rev_arrow=rxn_rev_arrow, reversible_arrow=rxn_reversible_arrow)
        gpr = extract(row, rxn_gpr_key)
        if "," in gpr or "+" in gpr:
            # break on ',' (which is or) then '+' (which is and)
            ors = ["( %s  )" % o.replace("+", " and ") if "+" in o else o
                   for o in gpr.split(",")]
            gpr = " or ".join(ors)
        rxn.gene_reaction_rule = gpr
        if rxn_lb_key is not None:
            rxn.lower_bound = float(row[rxn_lb_key])
        if rxn_ub_key is not None:
            rxn.upper_bound = float(row[rxn_ub_key])

    # fix upper and lower bounds if they include infinity
    inf = float("inf")
    max_lower_bound = max(abs(i) for i in m.reactions.list_attr("lower_bound")
                          if abs(i) < inf)
    max_upper_bound = max(abs(i) for i in m.reactions.list_attr("upper_bound")
                          if abs(i) < inf)
    inf_replace = max(max_upper_bound, max_lower_bound, 1000.) * 100

    def clip_inf(value):
        if value == inf:
            return inf_replace
        elif value == -inf:
            return -inf_replace
        else:
            return value
    for reaction in m.reactions:
        reaction.lower_bound = clip_inf(reaction.lower_bound)
        reaction.upper_bound = clip_inf(reaction.upper_bound)

    return m
