""" stereo functionality for reaction objects
"""
from typing import List

import yaml

import automol.chi
import automol.geom
import automol.graph
from automol.reac._0core import (
    Reaction,
    from_forward_reverse,
    ts_graph,
    set_ts_graph,
    mapping,
    without_stereo,
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
    tsg = automol.graph.without_stereo(ts_graph(rxn))
    tsg0 = automol.graph.without_dummy_atoms(tsg)

    srxns = []
    for stsg0 in automol.graph.expand_stereo(tsg0, enant=enant):
        # To keep both sets of dummy atoms, copy the stereo parities over
        stsg = automol.graph.set_stereo_parities(
            tsg, automol.graph.stereo_parities(stsg0)
        )
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
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)
    rcts_gra = automol.graph.union_from_sequence(rct_gras)
    prds_gra = automol.graph.union_from_sequence(prd_gras)

    # 1. Align products with reactants
    rgra = automol.graph.relabel(rcts_gra, mapping(rxn, "R", "T"))
    pgra = automol.graph.relabel(prds_gra, mapping(rxn, "P", "T"))

    # 2. Expand all TSs
    stsgs = automol.graph.ts.expand_stereo_for_reaction(ts_graph(rxn), rgra, pgra)

    # 3. Copy them into reaction objects
    srxns = []
    for stsg in stsgs:
        srxn = set_ts_graph(rxn, stsg)
        srxns.append(srxn)

    return tuple(srxns)


def from_old_string(rxn_str, one_indexed=True, stereo=True):
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
    ftsg0 = automol.graph.from_old_yaml_dictionary(ftsg_dct, one_indexed=one_indexed)

    rtsg_dct = {
        "atoms": yaml_dct["backward TS atoms"],
        "bonds": yaml_dct["backward TS bonds"],
    }
    rtsg0 = automol.graph.from_old_yaml_dictionary(rtsg_dct, one_indexed=one_indexed)

    ftsg = automol.graph.without_stereo(ftsg0)
    rtsg = automol.graph.without_stereo(rtsg0)
    rxn = from_forward_reverse(cla, ftsg, rtsg, rcts_keys, prds_keys)

    if stereo:
        rcts_gra = automol.graph.ts.reagents_graph_without_stereo(
            ftsg0, keep_stereo=True
        )
        prds_gra = automol.graph.ts.reagents_graph_without_stereo(
            rtsg0, keep_stereo=True
        )
        rct_gras = [
            automol.graph.subgraph(rcts_gra, ks, stereo=True) for ks in rcts_keys
        ]
        prd_gras = [
            automol.graph.subgraph(prds_gra, ks, stereo=True) for ks in prds_keys
        ]
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
    srxn = set_ts_graph(srxn, automol.raph.reflect(stsg))
    return srxn
