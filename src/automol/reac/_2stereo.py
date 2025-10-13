""" stereo functionality for reaction objects
"""

from collections.abc import Sequence

import yaml

from .. import graph
from ._0core import (
    Reaction,
    from_forward_reverse,
    from_string,
    product_graphs,
    products_keys,
    reactant_graphs,
    reactants_keys,
    reagent_mappings,
    set_products_keys,
    set_reactants_keys,
    set_ts_graph,
    ts_graph,
)


def expand_stereo(
    rxn: Reaction,
    symeq: bool = False,
    enant: bool = True,
    strained: bool = False,
    rct_gras: Sequence[object] | None = None,
    prd_gras: Sequence[object] | None = None,
) -> tuple[Reaction, ...]:
    """Expand all possible stereo assignments for the reactants and products of
    this reaction. Only includes possibilities that are mutually consistent
    with each other.

    :param rxn: a reaction object
    :param symeq: Include symmetrically equivalent TS stereoisomers?
    :param enant: Include all TS enantiomers, or only canonical ones?
    :param strained: Include stereoisomers which are highly strained?
    :param rct_gras: Optionally, expand stereo to match these reactant graphs
    :param rct_gras: Optionally, expand stereo to match these product graphs
    :returns: a sequence of reaction objects with stereo assignments
    """

    def _combined_reagents_graph(gras: Sequence[object], prod: bool = False) -> object:
        """Form a combined reactants or products graph from individual components.

        The resulting graph will match the atom keys of the TS, but will not contain any
        dummy atoms. These are not needed by the stereoexpansion filter, which removes
        them for comparison.

        :param gras: The reactant or product graphs
        :param prod: Do this for the products?, defaults to False
        :return: The combined reactants or products graph
        """
        # If none of the graphs have stereo, don't try to match them
        if not any(map(graph.has_stereo, gras)):
            return None

        gras = [graph.standard_keys(graph.without_dummy_atoms(g)) for g in gras]
        maps = reagent_mappings(rxn, prod=prod, shift_keys=False, dummy=False)
        gras = [graph.relabel(g, m) for g, m in zip(gras, maps, strict=True)]
        return graph.union_from_sequence(gras, check=True)

    rgra = None if rct_gras is None else _combined_reagents_graph(rct_gras, prod=False)
    pgra = None if prd_gras is None else _combined_reagents_graph(prd_gras, prod=True)

    tsg = graph.without_stereo(ts_graph(rxn))
    stsgs = graph.expand_stereo(
        tsg, symeq=symeq, enant=enant, strained=strained, rcts_gra=rgra, prds_gra=pgra
    )
    srxns = tuple(set_ts_graph(rxn, stsg) for stsg in stsgs)
    return srxns


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


def from_old_string(
    rxn_str: str,
    one_indexed: bool = True,
    stereo: bool = False,
    enant: bool = True,
    strained: bool = True,
) -> Reaction:
    """Write a reaction object to a string.

    :param rxn_str: string containing the reaction object
    :param one_indexed: parameter to store keys in one-indexing
    :param enant: If expanding stereo, include enantiomers?
    :param strained: If expanding stereo, include strained stereoisomers?
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
        rxns = expand_stereo(
            rxn,
            enant=enant,
            strained=strained,
            rct_gras=rct_gras,
            prd_gras=prd_gras,
        )
        rxn, *_ = rxns

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
