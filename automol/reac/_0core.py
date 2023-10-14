""" reaction classifiers and reaction-class-specific functions

Function arguments:
    Each function takes a list of reactant graphs and a list of product graphs.
    Note that the reactant graphs *cannot* have overlapping atom keys, and
    likewise for the product graphs. Otherwise, there would be no way to
    express the bonds broken and formed between reactants.
"""
import copy
import dataclasses
import itertools
from typing import Dict, List

import yaml

from automol import geom, graph, par, zmat
from automol.graph import ts
from automol.util import ZmatConv, zmat_conv


@dataclasses.dataclass
class Reaction:
    """Encodes essential information about a reaction, how the TS relates to
    the reactants and products. Also stores and allows handling of structural
    information in Cartesian geometry or z-matrix formats.

    :param ts_graph: The TS graph
    :type ts_graph: automol graph data structure
    :param reactants_keys: Keys to the atoms of each reactant in `ts_graph`
    :type reactants_keys: tuple[tuple[int]]
    :param reactants_keys: Keys to the atoms of each product in `ts_graph`
    :type products_keys: tuple[tuple[int]]
    :param class_: The reaction class
    :type class_: str
    :param structure_type: The structural information type ('zmat' or 'geom')
    :type structure_type: str
    :param ts_structure: The TS stucture, with keys matching the TS graph
    :type ts_structure: automol geom or zmat data structure
    :param reactant_structures: The reactant stuctures, with keys matching reactants
    :type reactant_structures: List[automol geom or zmat data structure]
    :param product_structures: The product stuctures, with keys matching products
    :type product_structures: List[automol geom or zmat data structure]
    :param ts_conversion_info: Z-matrix conversion info for the TS structure
    :type ts_conversion_info: ZmatConv
    :param reactants_conversion_info: Z-matrix conversion info for reactant structures
    :type reactants_conversion_info: ZmatConv
    :param products_conversion_info: Z-matrix conversion info for product structures
    :type products_conversion_info: ZmatConv
    """

    ts_graph: tuple
    reactants_keys: tuple
    products_keys: tuple
    class_: str
    structure_type: str = None
    ts_structure: tuple = None
    reactant_structures: tuple = None
    product_structures: tuple = None
    ts_conversion_info: ZmatConv = None
    reactants_conversion_info: List[ZmatConv] = None
    products_conversion_info: List[ZmatConv] = None

    def __repr__(self):
        """string representation of the object"""
        return string(self)


# # constructors
def from_forward_reverse(cla, ftsg, rtsg, rcts_keys, prds_keys) -> Reaction:
    """Construct a Reaction dataclass from forward and reverse TS graphs

    This function serves the (hopefully temporary) role of connecting the reaction
    finder to the new Reaction dataclass

    The forward TS graph matches the atom ordering of the products, while the
    reverse TS graph matches that of the reactants

    :param cla: The reaction class
    :type cla: str
    :param ftsg: The transition state graph, with reactant atom order
    :type ftsg: automol graph data structure
    :param rtsg: The reverse transition state graph, with product atom order
    :type rtsg: automol graph data structure
    :param rcts_keys: Keys to the reactant molecules in `ftsg`
    :type rcts_keys: tuple[tuple[int]]
    :param prds_keys: Keys to the product molecules in `rtsg`
    :type prds_keys: tuple[tuple[int]]
    :returns: A reaction object
    :rtype: Reaction
    """
    # Determine the reaction mapping
    rmap_dct = graph.isomorphism(ts.reverse(rtsg), ftsg, dummy=False, stereo=False)

    # Sort, so that the reagent orderings match the original input order
    rcts_keys = sorted(map(sorted, rcts_keys))
    prds_keys = sorted(map(sorted, prds_keys))

    # Reverse-map the products keys so they line up with the forward TS graph
    prds_keys = tuple(tuple(map(rmap_dct.__getitem__, ks)) for ks in prds_keys)

    return from_data(tsg=ftsg, rcts_keys=rcts_keys, prds_keys=prds_keys, cla=cla)


def from_data(
    tsg,
    rcts_keys,
    prds_keys,
    cla,
    ts_struc=None,
    rct_strucs=None,
    prd_strucs=None,
    ts_zc=None,
    rct_zcs=None,
    prd_zcs=None,
    struc_typ=None,
) -> Reaction:
    """Construct a Reaction dataclass from data

    :param tsg: The transition state graph
    :type tsg: automol graph data structure
    :param rcts_keys: Keys to the reactant molecules in `tsg`
    :type rcts_keys: tuple[tuple[int]]
    :param prds_keys: Keys to the product molecules in `tsg`
    :type prds_keys: tuple[tuple[int]]
    :param cla: The reaction class
    :type cla: str
    :param struc_type: The structural information type ('zmat' or 'geom')
    :type struc_type: str
    :param ts_struc: The TS stucture, with keys matching the TS graph
    :type ts_struc: automol geom or zmat data structure
    :param rct_strucs: The reactant stuctures, with keys matching reactants keys
    :type rct_strucs: List[automol geom or zmat data structure]
    :param prd_strucs: The product stuctures, with keys matching products keys
    :type prd_strucs: List[automol geom or zmat data structure]
    :param ts_zc: Z-matrix conversion info for the TS structure
    :type ts_zc: ZmatConv
    :param rct_zcs: Z-matrix conversion info for reactant structures
    :type rct_zcs: ZmatConv
    :param prd_zcs: Z-matrix conversion info for product structures
    :type prd_zcs: ZmatConv
    :returns: A reaction object
    :rtype: Reaction
    """

    # Check the TS graph...
    # If present, stereo information should be complete
    if graph.has_stereo(tsg):
        ste_keys = graph.stereogenic_keys(tsg)
        assert not ste_keys, f"TS graph has unassigned stereo at {ste_keys}:\n{tsg}"

    # Check the reactants and products keys
    rcts_keys = tuple(map(tuple, rcts_keys))
    prds_keys = tuple(map(tuple, prds_keys))

    # Check the reaction class
    assert par.is_reaction_class(cla) or cla is None, f"{cla} is not a reaction class"

    # Check the structures, if there are any
    struc_info = [ts_struc, rct_strucs, prd_strucs, struc_typ, ts_zc, rct_zcs, prd_zcs]
    if any(x is not None for x in struc_info):
        ttyp = _identify_sequence_structure_type([ts_struc])
        rtyp = _identify_sequence_structure_type(rct_strucs)
        ptyp = _identify_sequence_structure_type(prd_strucs)
        struc_typ = ttyp if struc_typ is None else struc_typ
        assert struc_typ is not None and ttyp == rtyp == ptyp == struc_typ, (
            f"Inconsistent structures:\ntype: {struc_typ}\nTS: {ts_struc}"
            f"\nreactants: {rct_strucs}\nproducts: {prd_strucs}"
        )

        round_ = geom.round_ if struc_typ == "geom" else zmat.round_
        symbs_ = geom.symbols if struc_typ == "geom" else zmat.symbols

        ts_struc = round_(ts_struc)
        rct_strucs = tuple(map(round_, rct_strucs))
        prd_strucs = tuple(map(round_, prd_strucs))

        # Check that TS graph and reagents keys match up with the respective structures
        s_ok = symbs_(ts_struc) == graph.symbols(tsg)
        assert s_ok, f"TS structure and graph don't match:\n{ts_struc}\n{tsg}"

        for struc, keys in zip(rct_strucs + prd_strucs, rcts_keys + prds_keys):
            s_ok = symbs_(struc) == symbs_(ts_struc, idxs=keys)
            assert s_ok, (
                f"Reagent structure and keys don't match TS graph...\n"
                f"Structure:\n{struc}\nKeys:\n{keys}\nTS graph:\n{tsg}"
            )

        # If present, check the z-matrix conversion info for correctness
        zc_info = [ts_zc, rct_zcs, prd_zcs]
        if any(x is not None for x in zc_info):
            assert all(x is not None for x in zc_info)
            assert len(rct_zcs) == len(rct_strucs)
            assert len(prd_zcs) == len(prd_strucs)

            rct_zcs = tuple(rct_zcs)
            prd_zcs = tuple(prd_zcs)

            # Check that the atom counts and dummy keys match each structure
            is_geom = struc_typ == "geom"
            count_ = geom.count if is_geom else zmat.count
            c_ok = zmat_conv.count(ts_zc, typ=struc_typ) == count_(ts_struc)
            d_ok = is_geom or zmat_conv.dummy_keys(ts_zc) == zmat.dummy_keys(ts_struc)
            assert (
                c_ok and d_ok
            ), f"TS conversion info doesn't match structure:{ts_zc}\n{ts_struc}"

            for zc_, struc in zip(rct_zcs + prd_zcs, rct_strucs + prd_strucs):
                c_ok = zmat_conv.count(zc_, typ=struc_typ) == count_(struc)
                d_ok = is_geom or zmat_conv.dummy_keys(zc_) == zmat.dummy_keys(struc)
                assert (
                    c_ok and d_ok
                ), f"Reagent conversion info doesn't match structure:\n{zc_}\n{struc}"

    return Reaction(
        ts_graph=tsg,
        reactants_keys=rcts_keys,
        products_keys=prds_keys,
        class_=cla,
        structure_type=struc_typ,
        ts_structure=ts_struc,
        reactant_structures=rct_strucs,
        product_structures=prd_strucs,
        ts_conversion_info=ts_zc,
        reactants_conversion_info=rct_zcs,
        products_conversion_info=prd_zcs,
    )


# # serialization
def from_string(rxn_str, one_indexed=True) -> Reaction:
    """Write a reaction object to a string

    TODO: Implement reading of structural info

    :param rxn_str: string containing the reaction object
    :type rxn_str: str
    :param one_indexed: parameter to store keys in one-indexing
    :type one_indexed: bool
    :returns: A reaction object
    :rtype: Reaction
    """
    yaml_dct = yaml.load(rxn_str, Loader=yaml.FullLoader)

    cla = yaml_dct["reaction class"]
    tsg = graph.from_yaml_data(yaml_dct, one_indexed=one_indexed)
    rcts_keys = tuple(map(tuple, yaml_dct["reactants keys"]))
    prds_keys = tuple(map(tuple, yaml_dct["products keys"]))

    if one_indexed:
        rcts_keys = [[k - 1 for k in ks] for ks in rcts_keys]
        prds_keys = [[k - 1 for k in ks] for ks in prds_keys]

    return from_data(cla=cla, tsg=tsg, rcts_keys=rcts_keys, prds_keys=prds_keys)


def string(rxn: Reaction, one_indexed=True) -> str:
    """Write a reaction object to a string

    :param rxn: the reaction object
    :type rxn: Reaction
    :param one_indexed: parameter to store keys in one-indexing
    :type one_indexed: bool
    :rtype: str
    """
    rcts_keys = list(map(list, reactants_keys(rxn)))
    prds_keys = list(map(list, products_keys(rxn)))

    if one_indexed:
        rcts_keys = [[k + 1 for k in ks] for ks in rcts_keys]
        prds_keys = [[k + 1 for k in ks] for ks in prds_keys]

    yaml_dct = graph.yaml_data(ts_graph(rxn), one_indexed=one_indexed)
    yaml_dct["reactants keys"] = list(map(list, rcts_keys))
    yaml_dct["products keys"] = list(map(list, prds_keys))
    yaml_dct["reaction class"] = class_(rxn)
    struc_typ = structure_type(rxn)
    if struc_typ is not None:
        yaml_ = geom.yaml_data if struc_typ == "geom" else zmat.yaml_data
        yaml_dct["structure type"] = struc_typ
        yaml_dct["TS structure"] = yaml_(ts_structure(rxn))
        yaml_dct["reactant structures"] = list(map(yaml_, reactant_structures(rxn)))
        yaml_dct["product structures"] = list(map(yaml_, product_structures(rxn)))

        ts_zc = ts_conversion_info(rxn)
        rct_zcs = reactants_conversion_info(rxn)
        prd_zcs = products_conversion_info(rxn)

        if ts_zc is not None:
            assert rct_zcs is not None and prd_zcs is not None
            yaml_dct["TS conversion info"] = zmat_conv.yaml_data(ts_zc)
            yaml_dct["reactants conversion info"] = list(
                map(zmat_conv.yaml_data, rct_zcs)
            )
            yaml_dct["products conversion info"] = list(
                map(zmat_conv.yaml_data, prd_zcs)
            )

    rxn_str = yaml.dump(yaml_dct, default_flow_style=None, sort_keys=False)
    return rxn_str


# # getters
def ts_graph(rxn: Reaction):
    """Get the TS graph of the reaction

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The TS graph
    :rtype: str
    """
    return rxn.ts_graph


def reactants_keys(rxn: Reaction) -> List[List[int]]:
    """Get the reactants keys of the reaction

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The reactants keys
    :rtype: List[List[int]]
    """
    return rxn.reactants_keys


def products_keys(rxn: Reaction) -> List[List[int]]:
    """Get the products keys of the reaction

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The products keys
    :rtype: List[List[int]]
    """
    return rxn.products_keys


def class_(rxn: Reaction) -> str:
    """Get the reaction class

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The reaction class
    :rtype: str
    """
    return rxn.class_


def structure_type(rxn: Reaction) -> str:
    """Get the type of the TS, reactant, and product structures (geom or zmat)

    Returns `None` if none have been set

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The structural information type ('zmat' or 'geom')
    :rtype: str
    """
    return rxn.structure_type


def ts_structure(rxn: Reaction):
    """Get the TS structure

    Returns `None` if none has been set

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The TS stucture, with keys matching the TS graph
    :rtype: automol geom or zmat data structure
    """
    return rxn.ts_structure


def reactant_structures(rxn: Reaction):
    """Get the reactant structures

    Returns `None` if none have been set

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The reactant stuctures, with keys matching reactants keys
    :rtype: List[automol geom or zmat data structure]
    """
    return rxn.reactant_structures


def product_structures(rxn: Reaction):
    """Get the product structures

    Returns `None` if none have been set

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The product stuctures, with keys matching products keys
    :rtype: List[automol geom or zmat data structure]
    """
    return rxn.product_structures


def ts_conversion_info(rxn: Reaction) -> ZmatConv:
    """Get z-matrix conversion info for the TS structure

    Returns `None` if none has been set

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: Z-matrix conversion info for the TS structure
    :rtype: ZmatConv
    """
    return rxn.ts_conversion_info


def reactants_conversion_info(rxn: Reaction) -> List[ZmatConv]:
    """Get z-matrix conversion info for the reactant structures

    Returns `None` if none have been set

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: Z-matrix conversion info for the reactant structures
    :rtype: List[ZmatConv]
    """
    return rxn.reactants_conversion_info


def products_conversion_info(rxn: Reaction) -> List[ZmatConv]:
    """Get z-matrix conversion info for the product structures

    Returns `None` if none have been set

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: Z-matrix conversion info for the product structures
    :rtype: List[ZmatConv]
    """
    return rxn.products_conversion_info


# # setters
def set_ts_graph(rxn: Reaction, tsg) -> Reaction:
    """Set the TS graph of the reaction

    :param rxn: The reaction object
    :type rxn: Reaction
    :param tsg: The TS graph
    :type tsg: automol graph data structure
    :returns: A new reaction object
    :rtype: Reaction
    """
    rxn = copy.deepcopy(rxn)
    rxn.ts_graph = tsg
    return rxn


def set_reactants_keys(rxn: Reaction, rcts_keys: List[List[int]]) -> Reaction:
    """Set the reactants keys of the reaction

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rcts_keys: The reactants keys
    :type rcts_keys: List[List[int]]
    :returns: A new reaction object
    :rtype: Reaction
    """
    rxn = copy.deepcopy(rxn)
    rxn.reactants_keys = rcts_keys
    return rxn


def set_products_keys(rxn: Reaction, prds_keys: List[List[int]]) -> Reaction:
    """Set the products keys of the reaction

    :param rxn: The reaction object
    :type rxn: Reaction
    :param prds_keys: The products keys
    :type prds_keys: List[List[int]]
    :returns: A new reaction object
    :rtype: Reaction
    """
    rxn = copy.deepcopy(rxn)
    rxn.products_keys = prds_keys
    return rxn


def set_reaction_class(rxn: Reaction, cla: str) -> Reaction:
    """Set the reaction class

    :param rxn: The reaction object
    :type rxn: Reaction
    :param cla: The reaction class
    :type cla: str
    :returns: A new reaction object
    :rtype: Reaction
    """
    rxn = copy.deepcopy(rxn)
    rxn.class_ = cla
    return rxn


def set_structures(
    rxn: Reaction,
    ts_struc,
    rct_strucs,
    prd_strucs,
    struc_typ=None,
    ts_zc: ZmatConv = None,
    rct_zcs: List[ZmatConv] = None,
    prd_zcs: List[ZmatConv] = None,
) -> Reaction:
    """Set the structures for the Reaction

    Structure type will be inferred if not passed in

    :param rxn: The reaction object
    :type rxn: Reaction
    :param ts_struc: The TS stuctures, with keys matching the TS graph
    :type ts_struc: automol geom or zmat data structure
    :param rct_strucs: The reactant stuctures, with keys matching reactants keys
    :type rct_strucs: List[automol geom or zmat data structure]
    :param prd_strucs: The product stuctures, with keys matching products keys
    :type prd_strucs: List[automol geom or zmat data structure]
    :param struc_typ: The structural information type ('zmat' or 'geom'),
        defaults to None
    :type struc_typ: str, optional
    :param ts_zc: Z-matrix conversion info for the TS structure
    :type ts_zc: ZmatConv
    :param rct_zcs: Z-matrix conversion info for reactant structures
    :type rct_zcs: ZmatConv
    :param prd_zcs: Z-matrix conversion info for product structures
    :type prd_zcs: ZmatConv
    :return: A new reaction object
    :rtype: Reaction
    """
    return from_data(
        tsg=ts_graph(rxn),
        rcts_keys=reactants_keys(rxn),
        prds_keys=products_keys(rxn),
        cla=class_(rxn),
        struc_typ=struc_typ,
        ts_struc=ts_struc,
        rct_strucs=rct_strucs,
        prd_strucs=prd_strucs,
        ts_zc=ts_zc,
        rct_zcs=rct_zcs,
        prd_zcs=prd_zcs,
    )


# # other
def relabel(rxn: Reaction, key_dct) -> Reaction:
    """Relabel keys in the TS graph

    :param rxn: the reaction object
    :type rxn: Reaction
    :param key_dct: A dictionary mapping current keys into new keys
    :type key_dct: dict[int: int]
    :returns: A relabeled reaction object
    :rtype: Reaction
    """
    tsg = graph.relabel(ts_graph(rxn), key_dct)
    rcts_keys = tuple(tuple(map(key_dct.__getitem__, ks)) for ks in reactants_keys(rxn))
    prds_keys = tuple(tuple(map(key_dct.__getitem__, ks)) for ks in products_keys(rxn))

    rxn = set_ts_graph(rxn, tsg)
    rxn = set_reactants_keys(rxn, rcts_keys)
    rxn = set_products_keys(rxn, prds_keys)
    return rxn


def reverse(rxn: Reaction) -> Reaction:
    """Obtains the reaction object for the reverse reaction

    TODO: Rename this reverse_without_structures and implement a general reversal
    function in _4struct.py that allows reversal with structures (including a flag for
    whether or not to keep the structures)

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: Reaction
    """
    return from_data(
        cla=par.reverse_reaction_class(class_(rxn)),
        tsg=graph.ts.reverse(ts_graph(rxn)),
        rcts_keys=products_keys(rxn),
        prds_keys=reactants_keys(rxn),
    )


def mapping(rxn: Reaction, inp: str, out: str) -> dict:
    """Determine a mapping from TS atoms to reactant atoms.

    Code:
        'T': TS keys
        'R': reactant keys
        'P': product keys

    :param rxn: A reaction object
    :type rxn: Reaction
    :param inp: Keys going into the mapping, 'T', 'R', or 'P'
    :type inp: str
    :param out: Keys coming out of the mapping, 'T', 'R', or 'P'
    :type out: str
    :returns: A dictionary mapping `inp` keys into `out` keys
    :rtype: dict
    """
    assert inp in ("T", "R", "P") and out in ("T", "R", "P")

    keys = sorted(graph.atom_keys(ts_graph(rxn)))
    rcts_keys = reactants_keys(rxn)
    prds_keys = products_keys(rxn)
    rct_keys = list(itertools.chain(*rcts_keys))
    prd_keys = list(itertools.chain(*prds_keys))

    key_dct = {
        "T": keys,
        "R": [rct_keys.index(k) if k in rct_keys else None for k in keys],
        "P": [prd_keys.index(k) if k in prd_keys else None for k in keys],
    }
    inp_keys = key_dct[inp]
    out_keys = key_dct[out]
    map_dct = dict(zip(inp_keys, out_keys))

    if None in map_dct:
        map_dct.pop(None)

    return map_dct


def reactant_mappings(
    rxn: Reaction, rev: bool = False, shift_keys: bool = False
) -> List[Dict[int, int]]:
    """Mappings from each reactant onto the TS

    :param rxn: A reaction object
    :type rxn: Reaction
    :param rev: Reverse the mapping? defaults to False
    :type rev: bool, optional
    :param shift_keys: Shift keys after first reagent, to prevent overlap? default False
    :type shift_keys: bool, optional
    :returns: The list of mappings for each reactant in order
    :rtype: List[Dict[int, int]]
    """
    start = 0
    map_dcts = []
    for keys in reactants_keys(rxn):
        map_dct = {k + start: v for k, v in enumerate(keys)}

        if rev:
            map_dct = dict(map(reversed, map_dct.items()))

        if shift_keys:
            start += len(keys)

        map_dcts.append(map_dct)

    return tuple(map_dcts)


def product_mappings(
    rxn: Reaction, rev: bool = False, shift_keys: bool = False
) -> List[Dict[int, int]]:
    """Mappings from each product onto the TS

    :param rxn: A reaction object
    :type rxn: Reaction
    :param rev: Reverse the mapping? defaults to False
    :type rev: bool, optional
    :param shift_keys: Shift keys after first reagent, to prevent overlap? default False
    :type shift_keys: bool, optional
    :returns: The list of mappings for each product in order
    :rtype: List[Dict[int, int]]
    """
    return reactant_mappings(reverse(rxn), rev=rev, shift_keys=shift_keys)


def reactant_graphs(rxn: Reaction, shift_keys: bool = False):
    """Obtain graphs of the reactants in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param shift_keys: Shift keys after first reagent, to prevent overlap? default False
    :type shift_keys: bool, optional
    :rtype: tuple of automol graph data structures
    """
    rcts_gra = ts.reactants_graph(ts_graph(rxn))
    rcts_keys = reactants_keys(rxn)

    # Extract subgraphs (TS keys)
    rct_gras = [graph.subgraph(rcts_gra, ks, stereo=True) for ks in rcts_keys]

    # Get mappings onto shifted or unshifted reactant keys
    map_dcts = reactant_mappings(rxn, rev=True, shift_keys=shift_keys)
    rct_gras = [graph.relabel(g, m) for g, m in zip(rct_gras, map_dcts)]
    return tuple(rct_gras)


def product_graphs(rxn: Reaction, shift_keys=False):
    """Obtain graphs of the products in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param shift_keys: Shift keys after first reagent, to prevent overlap? default False
    :type shift_keys: bool, optional
    :rtype: tuple of automol graph data structures
    """
    return reactant_graphs(reverse(rxn), shift_keys=shift_keys)


def reactants_graph(rxn: Reaction, key_order="R"):
    """Obtain a (single) graph of the reactants in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param key_order: The key order to use, 'T', 'R', or 'P'
    :type key_order: str, optional
    :rtype: automol graph data structure
    """
    map_dct = mapping(rxn, "T", key_order)
    rcts_gra = ts.reactants_graph(ts_graph(rxn))
    return graph.relabel(rcts_gra, map_dct, check=True)


def products_graph(rxn: Reaction, key_order="P"):
    """Obtain a (single) graph of the products in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param key_order: The key order to use, 'T', 'R', or 'P'
    :type key_order: str, optional
    :rtype: automol graph data structure
    """
    map_dct = mapping(rxn, "T", key_order)
    rcts_gra = ts.products_graph(ts_graph(rxn))
    return graph.relabel(rcts_gra, map_dct, check=True)


def standard_keys(rxn: Reaction) -> Reaction:
    """Standardize keys for the reaction object

    Ensures that keys follow zero-indexing in order with no skips

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: A copy of the reaction object, with standardized keys
    :rtype: Reaction
    """
    rct_keys = list(map(sorted, reactants_keys(rxn)))
    rct_key_dct = {k: i for i, k in enumerate(itertools.chain(*rct_keys))}

    rxn = relabel(rxn, rct_key_dct)
    return rxn


def without_stereo(rxn: Reaction) -> Reaction:
    """Remove all stereo information from the reaction object

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the reaction object, without stereo information
    :rtype: Reaction
    """
    tsg = ts_graph(rxn)
    rxn = set_ts_graph(rxn, graph.without_stereo(tsg))
    return rxn


def apply_zmatrix_conversion(rxn: Reaction, zc_: ZmatConv) -> Reaction:
    """Apply a z-matrix conversion (dummy insertion + reordering) to the reaction

    DEPRECATED -- this is taken care of by structure conversion

    This can be used to match a z-matrix after geometry -> z-matrix conversion

    :param rxn: A reaction object
    :type rxn: Reaction
    :param zc_: A z-matrix conversion
    :type zc_: ZmatConv
    :returns: The converted reaction object
    :rtype: Reaction
    """
    # 1. Transform the TS graph
    tsg = graph.apply_zmatrix_conversion(ts_graph(rxn), zc_)
    # 2. Relabel the keys
    rel_dct = zmat_conv.relabel_dict(zc_)
    rcts_keys = list(list(map(rel_dct.__getitem__, ks)) for ks in reactants_keys(rxn))
    prds_keys = list(list(map(rel_dct.__getitem__, ks)) for ks in products_keys(rxn))
    # 3. Insert the dummy keys in the right places
    # Temporary solution
    # Ultimately, I think we will need to have separate z-matrix conversions for each
    # reactant and product if we want everything to properly line up
    ins_dct = zmat_conv.insert_dict(zc_)
    for dummy_key, parent_key in sorted(ins_dct.items()):
        rkeys = next(ks for ks in rcts_keys if parent_key in ks)
        pkeys = next(ks for ks in prds_keys if parent_key in ks)
        rkeys.append(dummy_key)
        pkeys.append(dummy_key)

    rxn = set_ts_graph(rxn, tsg)
    rxn = set_reactants_keys(rxn, rcts_keys)
    rxn = set_products_keys(rxn, prds_keys)
    return rxn


def reverse_zmatrix_conversion(rxn: Reaction, zc_: ZmatConv) -> Reaction:
    """Reverse a z-matrix conversion (dummy insertion + reordering) to the reaction

    DEPRECATED -- this is taken care of by structure conversion

    This can be used to match the original geometry after conversio to z-matrix

    :param rxn: A converted reaction object
    :type rxn: Reaction
    :param zc_: A z-matrix conversion
    :type zc_: ZmatConv
    :returns: The original reaction object
    :rtype: Reaction
    """
    rxn = without_dummy_atoms(rxn)
    rel_dct = zmat_conv.relabel_dict(zc_, rev=True)
    rxn = relabel(rxn, rel_dct)
    return rxn


def without_dummy_atoms(rxn: Reaction) -> Reaction:
    """remove dummy atoms from the reactants or products"""
    tsg = ts_graph(rxn)
    dummy_keys = graph.atom_keys(tsg, symb="X")

    tsg = graph.without_dummy_atoms(tsg)
    rcts_keys = [[k for k in ks if k not in dummy_keys] for ks in reactants_keys(rxn)]
    prds_keys = [[k for k in ks if k not in dummy_keys] for ks in products_keys(rxn)]
    rxn = set_ts_graph(rxn, tsg)
    rxn = set_reactants_keys(rxn, rcts_keys)
    rxn = set_products_keys(rxn, prds_keys)
    return rxn


def relabel_for_geometry(rxn: Reaction) -> Reaction:
    """relabel the reaction object to correspond with a geometry converted
    from a z-matrix

    DEPRECATED

    :param rxn: the reaction object
    :param product: do this do the products instead of the reactants?
    """
    rxn = without_dummy_atoms(rxn)

    tsg = ts_graph(rxn)
    keys = sorted(graph.atom_keys(tsg))
    key_dct = dict(map(reversed, enumerate(keys)))
    rxn = relabel(rxn, key_dct)
    return rxn


def unique(rxns: List[Reaction]) -> List[Reaction]:
    """Get reactions with distinct TSs from a list with redundancies

    :param rxns: a sequence of reaction objects
    :returns: unique reaction objects
    """
    all_rxns = rxns
    rxns = []

    def isomorphic_(rxn1, rxn2):
        tsg1 = ts_graph(rxn1)
        tsg2 = ts_graph(rxn2)
        return graph.isomorphic(tsg1, tsg2, stereo=True)

    for rxn in all_rxns:
        if not any(isomorphic_(rxn, r) for r in rxns):
            rxns.append(rxn)

    return tuple(rxns)


def is_radical_radical(zrxn: Reaction) -> bool:
    """Is this a radical-radical reaction

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: boolean
    """
    _is_rad_rad = False
    rct_gras = reactant_graphs(zrxn)
    if len(rct_gras) == 2:
        rct_i, rct_j = rct_gras
        if graph.is_radical_species(rct_i) and graph.is_radical_species(rct_j):
            _is_rad_rad = True
    return _is_rad_rad


def filter_viable_reactions(rxns: List[Reaction]) -> List[Reaction]:
    """Filter a list of reactions to only include the viable ones

    Currently, filters reactions where:
     - One or more products has separated radical sites
     - The high-spin multiplicity of the products as a whole is greater than 3.

    :param rxns: a sequence of reactions
    :type rxns: tuple[Reaction]
    :returns: reactions with viable products
    :rtype: tuple[Reaction]
    """
    all_rxns = rxns
    rxns = []

    def _produces_separated_radical_sites(rxn):
        prd_gras = product_graphs(rxn)
        sep_rad = any(
            graph.has_separated_radical_sites(prd_gra) for prd_gra in prd_gras
        )
        return sep_rad

    def _high_spin_products(rxn):
        prd_gras = product_graphs(rxn)
        mult = sum(
            map(
                graph.maximum_spin_multiplicity,
                map(graph.kekule, prd_gras),
            )
        )
        # 4 allows for singlet+triplet or doublet+doublet products
        # 5 maybe best, allow for triplet+doublet (high-spin alkylrad+O2 HAbs)
        return mult > 4

    for rxn in all_rxns:
        # Check for separated radical sites
        sep_rad = _produces_separated_radical_sites(rxn)
        hi_spin = _high_spin_products(rxn)

        # Add more conditions here, as needed ...

        if not (sep_rad or hi_spin):
            rxns.append(rxn)

    return tuple(rxns)


# helpers
def _identify_sequence_structure_type(strucs):
    """Identify the type of a sequence of structures

    :param strucs: The structures (geom or zmat)
    :type strucs: List[automol geom or zmat data structure]
    """
    return (
        "geom"
        if all(map(geom.is_valid, strucs))
        else "zmat"
        if all(map(zmat.is_valid, strucs))
        else None
    )
