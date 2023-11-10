""" stereo functionality for reaction objects
"""
from typing import List

import yaml

from automol import graph
from automol.reac._0core import (
    Reaction,
    from_forward_reverse,
    mapping,
    set_ts_graph,
    ts_graph,
    without_stereo,
    from_string,
    reactant_graphs,
    product_graphs,
    reactants_keys,
    products_keys,
    set_reactants_keys,
    set_products_keys,
)


def expand_stereo(rxn: Reaction, enant=True) -> List[Reaction]:
    """Expand all possible stereo assignments for the reactants and products of
    this reaction. Only includes possibilities that are mutually consistent
    with each other.

    :param rxn: a reaction object
    :type rxn: Reaction
    :param enant: Include all enantiomers? Otherwise, includes only canonical
        enantiomer species and reactions.
    :type enant: bool
    :returns: a sequence of reaction objects with stereo assignments
    :rtype: List[Reaction]
    """
    tsg = graph.without_stereo(ts_graph(rxn))
    tsg0 = graph.without_dummy_atoms(tsg)

    srxns = []
    for stsg0 in graph.expand_stereo(tsg0, enant=enant):
        # To keep both sets of dummy atoms, copy the stereo parities over
        stsg = graph.set_stereo_parities(tsg, graph.stereo_parities(stsg0))
        srxn = set_ts_graph(rxn, stsg)

        srxns.append(srxn)

    srxns = tuple(srxns)
    return srxns


def expand_stereo_for_reaction(rxn: Reaction, rct_gras, prd_gras):
    """Expand stereo to be consistent with the reactants and products

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rcts_gra: A graph of the reactants
    :type rcts_gra: automol graph data structure
    :param prds_gra: A graph of the products
    :type prds_gra: automol graph data structure
    """
    rxn = without_stereo(rxn)

    # 0. Build combined reactants/products graphs
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)
    rcts_gra = graph.union_from_sequence(rct_gras)
    prds_gra = graph.union_from_sequence(prd_gras)

    # 1. Align products with reactants
    rgra = graph.relabel(rcts_gra, mapping(rxn, "R", "T"))
    pgra = graph.relabel(prds_gra, mapping(rxn, "P", "T"))

    # 2. Expand all TSs
    stsgs = graph.ts.expand_ts_stereo_for_reaction(ts_graph(rxn), rgra, pgra)

    # 3. Copy them into reaction objects
    srxns = []
    for stsg in stsgs:
        srxn = set_ts_graph(rxn, stsg)
        srxns.append(srxn)

    return tuple(srxns)


def from_string_transitional(rxn_str):
    """A transitional string reader which reads both old and new reaction strings

    :param rxn_str: string containing the (old or new) reaction object
    :type rxn_str: str
    """
    yaml_dct = yaml.load(rxn_str, Loader=yaml.FullLoader)
    if "forward TS atoms" in yaml_dct:
        rxn = from_old_string(rxn_str)
    else:
        assert "atoms" in yaml_dct, f"Reaction string not recognized:\n{rxn_str}"
        rxn = from_string(rxn_str)
    return rxn


def from_old_string(rxn_str, one_indexed=True, stereo=False):
    """Write a reaction object to a string

    :param rxn_str: string containing the reaction object
    :type rxn_str: str
    :param one_indexed: parameter to store keys in one-indexing
    :type one_indexed: bool
    :rtype: Reaction
    """
    yaml_dct = yaml.load(rxn_str, Loader=yaml.FullLoader)

    cla = yaml_dct["reaction class"]
    rcts_keys = yaml_dct["reactants keys"]
    prds_keys = yaml_dct["products keys"]

    if one_indexed:
        rcts_keys = [[k - 1 for k in ks] for ks in rcts_keys]
        prds_keys = [[k - 1 for k in ks] for ks in prds_keys]

    ftsg_dct = {
        "atoms": yaml_dct["forward TS atoms"],
        "bonds": yaml_dct["forward TS bonds"],
    }
    ftsg0 = graph.from_old_yaml_data(ftsg_dct, one_indexed=one_indexed)

    rtsg_dct = {
        "atoms": yaml_dct["backward TS atoms"],
        "bonds": yaml_dct["backward TS bonds"],
    }
    rtsg0 = graph.from_old_yaml_data(rtsg_dct, one_indexed=one_indexed)

    ftsg = graph.without_stereo(ftsg0)
    rtsg = graph.without_stereo(rtsg0)
    rxn = from_forward_reverse(cla, ftsg, rtsg, rcts_keys, prds_keys)

    # Handle dummy atoms
    if graph.has_dummy_atoms(ftsg):
        rcts_keys0 = reactants_keys(rxn)
        prds_keys0 = products_keys(rxn)
        rct_gras = reactant_graphs(rxn, key_order="T")
        prd_gras = product_graphs(rxn, key_order="T")

        rcts_keys = []
        for rct_keys0, rct_gra in zip(rcts_keys0, rct_gras):
            bad_dum_keys = graph.unneeded_dummy_atom_keys(rct_gra)
            rct_keys = [k for k in rct_keys0 if k not in bad_dum_keys]
            rcts_keys.append(rct_keys)

        prds_keys = []
        for prd_keys0, prd_gra in zip(prds_keys0, prd_gras):
            bad_dum_keys = graph.unneeded_dummy_atom_keys(prd_gra)
            prd_keys = [k for k in prd_keys0 if k not in bad_dum_keys]
            prds_keys.append(prd_keys)

        rxn = set_reactants_keys(rxn, rcts_keys)
        rxn = set_products_keys(rxn, prds_keys)

    if stereo:
        rcts_gra = graph.ts.reagents_graph_without_stereo(
            ftsg0, keep_stereo=True, dummy=True
        )
        prds_gra = graph.ts.reagents_graph_without_stereo(
            rtsg0, keep_stereo=True, dummy=True
        )
        rct_gras = [graph.subgraph(rcts_gra, ks, stereo=True) for ks in rcts_keys]
        prd_gras = [graph.subgraph(prds_gra, ks, stereo=True) for ks in prds_keys]
        # Now, work out the combined stereochemistry
        rxn = expand_stereo_for_reaction(rxn, rct_gras, prd_gras)[0]

    return rxn


def reflect(srxn: Reaction):
    """Reflect all graphs in this reaction, to obtain their mirror images

    :param srxn: a reaction object with stereo assignments
    :type srxn: Reaction
    :returns: a reflected reaction object
    """
    stsg = ts_graph(srxn)
    srxn = set_ts_graph(srxn, graph.reflect(stsg))
    return srxn
