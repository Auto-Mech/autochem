"""Reaction classifiers and reaction-class-specific functions.

Function arguments:
    Each function takes a list of reactant graphs and a list of product graphs.
    Note that the reactant graphs *cannot* have overlapping atom keys, and
    likewise for the product graphs. Otherwise, there would be no way to
    express the bonds broken and formed between reactants.
"""

import copy
import dataclasses
import itertools

import yaml
from phydat import ptab

from .. import const, geom, graph, zmat
from ..graph import ts
from ..util import ZmatConv, dict_, zmat_conv


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
    :type reactant_structures: list[automol geom or zmat data structure]
    :param product_structures: The product stuctures, with keys matching products
    :type product_structures: list[automol geom or zmat data structure]
    :param ts_conversion_info: Z-matrix conversion info for the TS structure
    :type ts_conversion_info: ZmatConv
    :param reactants_conversion_info: Z-matrix conversion info for reactant structures
    :type reactants_conversion_info: ZmatConv
    :param products_conversion_info: Z-matrix conversion info for product structures
    :type products_conversion_info: ZmatConv
    """

    ts_graph: tuple
    reactants_keys: tuple[tuple[int, ...], ...]
    products_keys: tuple[tuple[int, ...], ...]
    class_: str
    structure_type: str = None
    ts_structure: tuple = None
    reactant_structures: tuple = None
    product_structures: tuple = None
    ts_conversion_info: ZmatConv = None
    reactants_conversion_info: list[ZmatConv] = None
    products_conversion_info: list[ZmatConv] = None

    def __repr__(self):
        """string representation of the object"""
        return string(self)


# # constructors
def from_data(
    tsg,
    rcts_keys,
    prds_keys,
    cla=const.ReactionClass.UNCLASSIFIED,
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
    :type rct_strucs: list[automol geom or zmat data structure]
    :param prd_strucs: The product stuctures, with keys matching products keys
    :type prd_strucs: list[automol geom or zmat data structure]
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
        ste_keys = graph.unassigned_stereocenter_keys(tsg)
        assert not ste_keys, f"TS graph has unassigned stereo at {ste_keys}:\n{tsg}"

    # Check the reactants and products keys
    rcts_keys = tuple(map(tuple, rcts_keys))
    prds_keys = tuple(map(tuple, prds_keys))

    # Check the reaction class
    assert (
        const.ReactionClass.is_defined(cla) or cla is None
    ), f"{cla} is not a reaction class"

    # Check the structures, if there are any
    struc_info = [ts_struc, rct_strucs, prd_strucs]
    if any(x is not None for x in struc_info):
        assert all(x is not None for x in struc_info)
        assert len(rct_strucs) == len(rcts_keys)
        assert len(prd_strucs) == len(prds_keys)

        rct_strucs = tuple(rct_strucs)
        prd_strucs = tuple(prd_strucs)

        # Infer the structural type
        ttyp = _identify_sequence_structure_type([ts_struc])
        rtyp = _identify_sequence_structure_type(rct_strucs)
        ptyp = _identify_sequence_structure_type(prd_strucs)
        struc_typ = ttyp if struc_typ is None else struc_typ
        assert struc_typ is not None and ttyp == rtyp == ptyp == struc_typ, (
            f"Inconsistent structures:\ntype: {struc_typ}\nTS: {ts_struc}"
            f"\nreactants: {rct_strucs}\nproducts: {prd_strucs}"
        )

        symbs_ = geom.symbols if struc_typ == "geom" else zmat.symbols

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
        assert len(rct_zcs) == len(rcts_keys)
        assert len(prd_zcs) == len(prds_keys)

        rct_zcs = tuple(rct_zcs)
        prd_zcs = tuple(prd_zcs)

        # Check that the atom counts and dummy keys match each structure
        no_zmats = struc_typ == "geom" or all(x is None for x in struc_info)
        c_ok = zmat_conv.count(ts_zc, typ=struc_typ) == graph.count(tsg)
        d_ok = no_zmats or zmat_conv.dummy_keys(ts_zc) == zmat.dummy_keys(ts_struc)
        assert (
            c_ok and d_ok
        ), f"TS conversion info doesn't match structure:{ts_zc}\n{ts_struc}"

        zcs = rct_zcs + prd_zcs
        keys_lst = rcts_keys + prds_keys
        zmas = [None] * len(zcs) if no_zmats else rct_strucs + prd_strucs
        for zc_, keys, zma in zip(zcs, keys_lst, zmas):
            c_ok = zmat_conv.count(zc_, typ=struc_typ) == len(keys)
            d_ok = no_zmats or zmat_conv.dummy_keys(zc_) == zmat.dummy_keys(zma)
            assert (
                c_ok and d_ok
            ), f"Reagent conversion info doesn't match structure:\n{zc_}\n{zma}"

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
    rmap_dct = graph.isomorphism(
        ts.reverse(rtsg), ftsg, dummy=False, stereo=False, unstable=False
    )

    # Sort, so that the reagent orderings match the original input order
    rcts_keys = sorted(map(sorted, rcts_keys))
    prds_keys = sorted(map(sorted, prds_keys))

    # Reverse-map the products keys so they line up with the forward TS graph
    prds_keys = tuple(tuple(map(rmap_dct.__getitem__, ks)) for ks in prds_keys)

    return from_data(tsg=ftsg, rcts_keys=rcts_keys, prds_keys=prds_keys, cla=cla)


# # serialization
def from_string(rxn_str) -> Reaction:
    """Write a reaction object to a string

    :param rxn_str: string containing the reaction object
    :type rxn_str: str
    :returns: A reaction object
    :rtype: Reaction
    """
    yaml_dct = yaml.load(rxn_str, Loader=yaml.FullLoader)

    tsg = graph.from_yaml_data(yaml_dct)
    rcts_keys = tuple(map(tuple, yaml_dct["reactants keys"]))
    prds_keys = tuple(map(tuple, yaml_dct["products keys"]))
    cla = const.ReactionClass(yaml_dct["reaction class"])

    struc_typ = yaml_dct.get("structure type")
    ts_struc = yaml_dct.get("ts structure")
    rct_strucs = yaml_dct.get("reactant structures")
    prd_strucs = yaml_dct.get("product structures")
    ts_zc = yaml_dct.get("ts conversion info")
    rct_zcs = yaml_dct.get("reactants conversion info")
    prd_zcs = yaml_dct.get("products conversion info")

    struc_info = [ts_struc, rct_strucs, prd_strucs]
    if any(x is not None for x in struc_info):
        assert struc_typ is not None and all(x is not None for x in struc_info)

        from_yaml_ = geom.from_yaml_data if struc_typ == "geom" else zmat.from_yaml_data
        ts_struc = from_yaml_(ts_struc)
        rct_strucs = list(map(from_yaml_, rct_strucs))
        prd_strucs = list(map(from_yaml_, prd_strucs))

    zc_info = [ts_zc, rct_zcs, prd_zcs]
    if any(x is not None for x in zc_info):
        assert all(x is not None for x in zc_info)

        ts_zc = zmat_conv.from_yaml_data(ts_zc)
        rct_zcs = list(map(zmat_conv.from_yaml_data, rct_zcs))
        prd_zcs = list(map(zmat_conv.from_yaml_data, prd_zcs))

    return from_data(
        cla=cla,
        tsg=tsg,
        rcts_keys=rcts_keys,
        prds_keys=prds_keys,
        ts_struc=ts_struc,
        rct_strucs=rct_strucs,
        prd_strucs=prd_strucs,
        ts_zc=ts_zc,
        rct_zcs=rct_zcs,
        prd_zcs=prd_zcs,
        struc_typ=struc_typ,
    )


def string(rxn: Reaction) -> str:
    """Write a reaction object to a string

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: str
    """
    rcts_keys = list(map(list, reactants_keys(rxn)))
    prds_keys = list(map(list, products_keys(rxn)))

    yaml_dct = graph.yaml_data(ts_graph(rxn))
    yaml_dct["reactants keys"] = list(map(list, rcts_keys))
    yaml_dct["products keys"] = list(map(list, prds_keys))
    yaml_dct["reaction class"] = str(class_(rxn))

    struc_typ = structure_type(rxn)
    ts_struc = ts_structure(rxn)
    rct_strucs = reactant_structures(rxn)
    prd_strucs = product_structures(rxn)
    ts_zc = ts_conversion_info(rxn)
    rct_zcs = reactants_conversion_info(rxn)
    prd_zcs = products_conversion_info(rxn)

    struc_info = [ts_struc, rct_strucs, prd_strucs]
    if any(x is not None for x in struc_info):
        assert struc_typ is not None and all(x is not None for x in struc_info)

        yaml_ = geom.yaml_data if struc_typ == "geom" else zmat.yaml_data
        ts_struc = yaml_(ts_struc)
        rct_strucs = list(map(yaml_, rct_strucs))
        prd_strucs = list(map(yaml_, prd_strucs))

    zc_info = [ts_zc, rct_zcs, prd_zcs]
    if any(x is not None for x in zc_info):
        assert all(x is not None for x in zc_info)

        ts_zc = zmat_conv.yaml_data(ts_zc)
        rct_zcs = list(map(zmat_conv.yaml_data, rct_zcs))
        prd_zcs = list(map(zmat_conv.yaml_data, prd_zcs))

    yaml_dct["structure type"] = struc_typ
    yaml_dct["ts structure"] = ts_struc
    yaml_dct["reactant structures"] = rct_strucs
    yaml_dct["product structures"] = prd_strucs
    yaml_dct["ts conversion info"] = ts_zc
    yaml_dct["reactants conversion info"] = rct_zcs
    yaml_dct["products conversion info"] = prd_zcs

    rxn_str = yaml.dump(yaml_dct, default_flow_style=None, sort_keys=False)
    return rxn_str


# # getters
def ts_graph(rxn: Reaction):
    """Get the TS graph of the reaction.

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The TS graph
    :rtype: str
    """
    return rxn.ts_graph


def reactants_keys(rxn: Reaction, dummy: bool = True) -> tuple[tuple[int, ...], ...]:
    """Get the reactants keys of the reaction.

    :param rxn: The reaction object
    :param dummy: Include dummy atoms?
    :returns: The reactants keys
    """
    rcts_keys = rxn.reactants_keys
    if not dummy:
        sym_dct = graph.atom_symbols(ts_graph(rxn))
        rcts_keys = [
            [k for k in ks if ptab.to_number(sym_dct.get(k))] for ks in rcts_keys
        ]
    return tuple(map(tuple, rcts_keys))


def products_keys(rxn: Reaction, dummy: bool = True) -> tuple[tuple[int, ...], ...]:
    """Get the products keys of the reaction.

    :param rxn: The reaction object
    :param dummy: Include dummy atoms?
    :returns: The products keys
    """
    prds_keys = rxn.products_keys
    if not dummy:
        sym_dct = graph.atom_symbols(ts_graph(rxn))
        prds_keys = [
            [k for k in ks if ptab.to_number(sym_dct.get(k))] for ks in prds_keys
        ]
    return tuple(map(tuple, prds_keys))


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
    :rtype: list[automol geom or zmat data structure]
    """
    return rxn.reactant_structures


def product_structures(rxn: Reaction):
    """Get the product structures

    Returns `None` if none have been set

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The product stuctures, with keys matching products keys
    :rtype: list[automol geom or zmat data structure]
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


def reactants_conversion_info(rxn: Reaction) -> list[ZmatConv]:
    """Get z-matrix conversion info for the reactant structures

    Returns `None` if none have been set

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: Z-matrix conversion info for the reactant structures
    :rtype: list[ZmatConv]
    """
    return rxn.reactants_conversion_info


def products_conversion_info(rxn: Reaction) -> list[ZmatConv]:
    """Get z-matrix conversion info for the product structures

    Returns `None` if none have been set

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: Z-matrix conversion info for the product structures
    :rtype: list[ZmatConv]
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
    return from_data(
        tsg=tsg,
        rcts_keys=reactants_keys(rxn),
        prds_keys=products_keys(rxn),
        cla=class_(rxn),
        ts_struc=ts_structure(rxn),
        rct_strucs=reactant_structures(rxn),
        prd_strucs=product_structures(rxn),
        struc_typ=structure_type(rxn),
    )


def set_reactants_keys(rxn: Reaction, rcts_keys: list[list[int]]) -> Reaction:
    """Set the reactants keys of the reaction

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rcts_keys: The reactants keys
    :type rcts_keys: list[list[int]]
    :returns: A new reaction object
    :rtype: Reaction
    """
    return from_data(
        tsg=ts_graph(rxn),
        rcts_keys=rcts_keys,
        prds_keys=products_keys(rxn),
        cla=class_(rxn),
        ts_struc=ts_structure(rxn),
        rct_strucs=reactant_structures(rxn),
        prd_strucs=product_structures(rxn),
        struc_typ=structure_type(rxn),
    )


def set_products_keys(rxn: Reaction, prds_keys: list[list[int]]) -> Reaction:
    """Set the products keys of the reaction

    :param rxn: The reaction object
    :type rxn: Reaction
    :param prds_keys: The products keys
    :type prds_keys: list[list[int]]
    :returns: A new reaction object
    :rtype: Reaction
    """
    return from_data(
        tsg=ts_graph(rxn),
        rcts_keys=reactants_keys(rxn),
        prds_keys=prds_keys,
        cla=class_(rxn),
        ts_struc=ts_structure(rxn),
        rct_strucs=reactant_structures(rxn),
        prd_strucs=product_structures(rxn),
        struc_typ=structure_type(rxn),
    )


def set_reaction_class(rxn: Reaction, cla: str) -> Reaction:
    """Set the reaction class.

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


def reset_conversion_info(rxn: Reaction) -> Reaction:
    """Reset z-matrix conversion info.

    For `geom` structures, this wipes out the conversion info.

    For `zmat` structures, this replaces the conversion info with simple dummy atom
    removal (no reordering).

    :param rxn: Reaction
    :return: Reaction
    """
    ts_zc = rct_zcs = prd_zcs = None

    if structure_type(rxn) == "zmat":
        ts_zc = zmat.conversion_info(ts_structure(rxn))
        rct_zcs = list(map(zmat.conversion_info, reactant_structures(rxn)))
        prd_zcs = list(map(zmat.conversion_info, product_structures(rxn)))

    return from_data(
        tsg=ts_graph(rxn),
        rcts_keys=reactants_keys(rxn),
        prds_keys=products_keys(rxn),
        cla=class_(rxn),
        struc_typ=structure_type(rxn),
        ts_struc=ts_structure(rxn),
        rct_strucs=reactant_structures(rxn),
        prd_strucs=product_structures(rxn),
        ts_zc=ts_zc,
        rct_zcs=rct_zcs,
        prd_zcs=prd_zcs,
    )


def set_structures(
    rxn: Reaction,
    ts_struc,
    rct_strucs,
    prd_strucs,
    struc_typ=None,
    ts_zc: ZmatConv = None,
    rct_zcs: list[ZmatConv] = None,
    prd_zcs: list[ZmatConv] = None,
) -> Reaction:
    """Set the structures for the Reaction.

    Structure type will be inferred if not passed in

    :param rxn: The reaction object
    :type rxn: Reaction
    :param ts_struc: The TS stuctures, with keys matching the TS graph
    :type ts_struc: automol geom or zmat data structure
    :param rct_strucs: The reactant stuctures, with keys matching reactants keys
    :type rct_strucs: list[automol geom or zmat data structure]
    :param prd_strucs: The product stuctures, with keys matching products keys
    :type prd_strucs: list[automol geom or zmat data structure]
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


def update_structures(
    rxn: Reaction,
    ts_struc=None,
    rct_strucs=None,
    prd_strucs=None,
    struc_typ=None,
    ts_zc: ZmatConv = None,
    rct_zcs: list[ZmatConv] = None,
    prd_zcs: list[ZmatConv] = None,
) -> Reaction:
    """Update structures or structural information for the Reaction

    :param rxn: The reaction object
    :type rxn: Reaction
    :param ts_struc: The TS stuctures, with keys matching the TS graph
    :type ts_struc: automol geom or zmat data structure
    :param rct_strucs: The reactant stuctures, with keys matching reactants keys
    :type rct_strucs: list[automol geom or zmat data structure]
    :param prd_strucs: The product stuctures, with keys matching products keys
    :type prd_strucs: list[automol geom or zmat data structure]
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
        struc_typ=structure_type(rxn) if struc_typ is None else struc_typ,
        ts_struc=ts_structure(rxn) if ts_struc is None else ts_struc,
        rct_strucs=reactant_structures(rxn) if rct_strucs is None else rct_strucs,
        prd_strucs=product_structures(rxn) if prd_strucs is None else prd_strucs,
        ts_zc=ts_conversion_info(rxn) if ts_zc is None else ts_zc,
        rct_zcs=reactants_conversion_info(rxn) if rct_zcs is None else rct_zcs,
        prd_zcs=products_conversion_info(rxn) if prd_zcs is None else prd_zcs,
    )


# # structures
def has_structures(rxn: Reaction, complete: bool = True) -> bool:
    """Does this reaction object have structures?

    :param rxn: A reaction object
    :type rxn: Reaction
    :param complete: Require all information to be present; defaults to True
    :type complete: bool, optional
    :return: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    strucs_lst = [ts_structure(rxn), reactant_structures(rxn), product_structures(rxn)]
    if complete:
        strucs_lst.append(structure_type(rxn))

    results = [x is not None for x in strucs_lst]
    return all(results) if complete else any(results)


def has_conversion_information(rxn: Reaction, complete: bool = True) -> bool:
    """Does this reaction object have structures?

    :param rxn: A reaction object
    :type rxn: Reaction
    :param complete: Require all information to be present; defaults to True
    :type complete: bool, optional
    :return: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    struc_info_lst = [
        ts_conversion_info(rxn),
        reactants_conversion_info(rxn),
        products_conversion_info(rxn),
    ]
    if complete:
        struc_info_lst.append(structure_type(rxn))

    results = [x is not None for x in struc_info_lst]
    return all(results) if complete else any(results)


def without_structures(rxn: Reaction, keep_info: bool = True) -> Reaction:
    """Remove structures from the reaction object

    :param rxn: A reaction object
    :type rxn: Reaction
    :param keep_zc: Keep structural conversion info?, defaults to True
    :type keep_zc: bool, optional
    :return: The reaction object, without structural information
    :rtype: Reaction
    """
    return from_data(
        tsg=ts_graph(rxn),
        rcts_keys=reactants_keys(rxn),
        prds_keys=products_keys(rxn),
        cla=class_(rxn),
        struc_typ=structure_type(rxn) if keep_info else None,
        ts_zc=ts_conversion_info(rxn) if keep_info else None,
        rct_zcs=reactants_conversion_info(rxn) if keep_info else None,
        prd_zcs=products_conversion_info(rxn) if keep_info else None,
    )


# # other
def relabel(rxn: Reaction, key_dct: dict[int, int], struc: bool = False) -> Reaction:
    """Relabel keys in the TS graph

    :param rxn: the reaction object
    :type rxn: Reaction
    :param key_dct: A dictionary mapping current keys into new keys
    :type key_dct: dict[int: int]
    :param struc: Keep the structures and relabel them?, defaults to False
    :type struc: bool, optional
    :returns: A relabeled reaction object
    :rtype: Reaction
    """
    tsg = graph.relabel(ts_graph(rxn), key_dct)
    rcts_keys = tuple(tuple(map(key_dct.get, ks)) for ks in reactants_keys(rxn))
    prds_keys = tuple(tuple(map(key_dct.get, ks)) for ks in products_keys(rxn))

    if not struc:
        rxn = without_structures(rxn)

    rxn = copy.copy(rxn)
    rxn.ts_graph = tsg
    rxn.reactants_keys = rcts_keys
    rxn.products_keys = prds_keys

    ts_struc = ts_structure(rxn)
    ts_zc = ts_conversion_info(rxn)
    if struc:
        struc_typ = structure_type(rxn)

        if struc_typ == "zmat":
            raise NotImplementedError("Cannot relabel z-matrix structures")

        if ts_struc is not None:
            ts_struc = geom.reorder(ts_struc, key_dct)

        if ts_zc is not None:
            ts_zc = zmat_conv.relabel(ts_zc, key_dct, "geom")

    return update_structures(rxn, ts_struc=ts_struc, ts_zc=ts_zc)


def reverse_without_recalculating(rxn: Reaction, struc: bool = True) -> Reaction:
    """Get the reaction object for the reverse reaction

    :param rxn: A reaction object
    :type rxn: Reaction
    :param struc: Keep structures?
    :type struc: bool, optional
    :returns: The reversed reaction object
    :rtype: Reaction
    """
    return from_data(
        cla=const.ReactionClass.reverse(class_(rxn)),
        tsg=graph.ts.reverse(ts_graph(rxn)),
        rcts_keys=products_keys(rxn),
        prds_keys=reactants_keys(rxn),
        ts_struc=ts_structure(rxn) if struc else None,
        rct_strucs=product_structures(rxn) if struc else None,
        prd_strucs=reactant_structures(rxn) if struc else None,
        ts_zc=ts_conversion_info(rxn) if struc else None,
        rct_zcs=products_conversion_info(rxn) if struc else None,
        prd_zcs=reactants_conversion_info(rxn) if struc else None,
        struc_typ=structure_type(rxn) if struc else None,
    )


def mapping(rxn: Reaction, inp: str, out: str, keep_none: bool = False) -> dict:
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
    :param keep_none: Keep entries in the mapping with none values? defaults to False
    :type keep_none: bool, optional
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
    map_dct = dict(zip(inp_keys, out_keys, strict=True))

    if None in map_dct:
        map_dct.pop(None)

    if not keep_none:
        map_dct = dict_.filter_by_value(map_dct, lambda x: x is not None)

    return map_dct


def reactant_mappings(
    rxn: Reaction, rev: bool = False, shift_keys: bool = False, dummy: bool = True
) -> tuple[dict[int, int], ...]:
    """Determine mappings from each reactant onto the TS.

    :param rxn: A reaction object
    :param rev: Reverse the mapping? defaults to False
    :param shift_keys: Shift keys after first reagent, to prevent overlap? default False
    :param dummy: Include dummy atoms?
    :returns: The list of mappings for each reactant in order
    """
    return reagent_mappings(
        rxn, prod=False, rev=rev, shift_keys=shift_keys, dummy=dummy
    )


def product_mappings(
    rxn: Reaction, rev: bool = False, shift_keys: bool = False, dummy: bool = True
) -> tuple[dict[int, int], ...]:
    """Determine mappings from each product onto the TS.

    :param rxn: A reaction object
    :param rev: Reverse the mapping? defaults to False
    :param shift_keys: Shift keys after first reagent, to prevent overlap? default False
    :param dummy: Include dummy atoms?
    :returns: The list of mappings for each product in order
    """
    return reagent_mappings(rxn, prod=True, rev=rev, shift_keys=shift_keys, dummy=dummy)


def reagent_mappings(
    rxn: Reaction,
    prod: bool = False,
    rev: bool = False,
    shift_keys: bool = False,
    dummy: bool = True,
) -> tuple[dict[int, int], ...]:
    """Determine mappings from each reactant or product onto the TS.

    :param rxn: A reaction object
    :param prod: Do this for the products?
    :param rev: Reverse the mapping? defaults to False
    :param shift_keys: Shift keys after first reagent, to prevent overlap? default False
    :param dummy: Include dummy atoms?
    :returns: The list of mappings for each reactant in order
    """
    if prod:
        rxn = reverse_without_recalculating(rxn)

    start = 0
    map_dcts = []
    for keys in reactants_keys(rxn, dummy=dummy):
        map_dct = {k + start: v for k, v in enumerate(keys)}

        if rev:
            map_dct = dict(map(reversed, map_dct.items()))

        if shift_keys:
            start += len(keys)

        map_dcts.append(map_dct)

    return tuple(map_dcts)


def reactant_graphs(
    rxn: Reaction, stereo: bool = True, key_order: str = "R", shift_keys: bool = False
):
    """Obtain graphs of the reactants in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param stereo: Include stereo? defaults to True
    :type stereo: bool, optional
    :param key_order: The key order to use, 'T' or 'R'
    :type key_order: str, optional
    :param shift_keys: Shift keys after first reagent, to prevent overlap? default False
        (Only has an effect on 'R' keys)
    :type shift_keys: bool, optional
    :rtype: tuple of automol graph data structures
    """
    assert key_order in ("T", "R")

    rcts_gra = ts.reactants_graph(ts_graph(rxn), stereo=stereo)
    rcts_keys = reactants_keys(rxn)

    # Extract subgraphs (TS keys)
    rct_gras = [graph.subgraph(rcts_gra, ks, stereo=True) for ks in rcts_keys]

    if key_order == "R":
        # Get mappings onto shifted or unshifted reactant keys
        map_dcts = reactant_mappings(rxn, rev=True, shift_keys=shift_keys)
        rct_gras = [graph.relabel(g, m) for g, m in zip(rct_gras, map_dcts)]
    return tuple(rct_gras)


def product_graphs(
    rxn: Reaction, stereo: bool = True, key_order: str = "P", shift_keys=False
):
    """Obtain graphs of the products in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param stereo: Include stereo? defaults to True
    :type stereo: bool, optional
    :param key_order: The key order to use, 'T' or 'P'
    :type key_order: str, optional
    :param shift_keys: Shift keys after first reagent, to prevent overlap? default False
        (Only has an effect on 'P' keys)
    :type shift_keys: bool, optional
    :type shift_keys: bool, optional
    :rtype: tuple of automol graph data structures
    """
    assert key_order in ("T", "P")

    rev_key_order = "T" if key_order == "T" else "R"

    return reactant_graphs(
        reverse_without_recalculating(rxn),
        stereo=stereo,
        key_order=rev_key_order,
        shift_keys=shift_keys,
    )


def reactants_graph(rxn: Reaction, stereo: bool = True, key_order="R"):
    """Obtain a (single) graph of the reactants in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param stereo: Include stereo? defaults to True
    :type stereo: bool, optional
    :param key_order: The key order to use, 'T', 'R', or 'P'
    :type key_order: str, optional
    :rtype: automol graph data structure
    """
    map_dct = mapping(rxn, "T", key_order)
    rcts_gra = ts.reactants_graph(ts_graph(rxn), stereo=stereo)
    return graph.relabel(rcts_gra, map_dct, check=True)


def products_graph(rxn: Reaction, stereo: bool = True, key_order="P"):
    """Obtain a (single) graph of the products in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param stereo: Include stereo? defaults to True
    :type stereo: bool, optional
    :param key_order: The key order to use, 'T', 'R', or 'P'
    :type key_order: str, optional
    :rtype: automol graph data structure
    """
    map_dct = mapping(rxn, "T", key_order)
    rcts_gra = ts.products_graph(ts_graph(rxn), stereo=stereo)
    return graph.relabel(rcts_gra, map_dct, check=True)


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


def apply_zmatrix_conversion(
    rxn: Reaction,
    ts_zc: ZmatConv = None,
    rct_zcs: list[ZmatConv] = None,
    prd_zcs: list[ZmatConv] = None,
) -> Reaction:
    """Apply a z-matrix conversion (dummy insertion + reordering) to the reaction

    This can be used to match a z-matrix after geometry -> z-matrix conversion

    :param rxn: A reaction object
    :type rxn: Reaction
    :param ts_zc: A z-matrix conversion for the TS graph, defaults to None
    :type ts_zc: ZmatConv
    :param rct_zcs: Z-matrix conversions for reactant structures, defaults to None
    :type rct_zcs: list[ZmatConv], optional
    :param prd_zcs: Z-matrix conversions for product structures, defaults to None
    :type prd_zcs: list[ZmatConv], optional
    :returns: The converted reaction object
    :rtype: Reaction
    """
    assert not has_structures(
        rxn, complete=False
    ), f"This should not be done for a reaction with structures:\n{rxn}"

    # Read in pre-existing z-matrix conversion info, if there is some
    ts_zc = ts_conversion_info(rxn) if ts_zc is None else ts_zc
    rct_zcs = reactants_conversion_info(rxn) if rct_zcs is None else rct_zcs
    prd_zcs = products_conversion_info(rxn) if prd_zcs is None else prd_zcs

    assert not any(
        x is None for x in [ts_zc, rct_zcs, prd_zcs]
    ), f"Missing conversion info:\nts_zc:{ts_zc}\nrct_zcs:{rct_zcs}\nprd_zcs:{prd_zcs}"

    # 0. Read in the current TS graph
    gtsg = ts_graph(rxn)
    assert not graph.has_dummy_atoms(
        gtsg
    ), f"Cannot apply z-matrix conversion to a graph with dummy atoms:{gtsg}"

    # 1. Transform the TS graph
    ztsg = graph.apply_zmatrix_conversion(gtsg, ts_zc)

    # 2. Transform the reagent keys
    def _get_reagent_zmat_keys(rgts_gkeys, rgt_zcs, rgt_gkey_dcts):
        """
        :param rgts_gkeys: Keys for each reagent realtive to the TS geometry
        :param rgt_zcs: Z-matrix conversions for each reagent
        :param rgt_gkey_dcts: Mappings from TS gkeys onto each reagent's gkeys
        """
        rgts_zkeys = []
        for rgt_gkeys, rgt_zc, rgt_gkey_dct in zip(rgts_gkeys, rgt_zcs, rgt_gkey_dcts):
            # 1. {TS zkeys: TS gkeys} (extract the subset matching reagent gkeys)
            zc_ = zmat_conv.subset(ts_zc, rgt_gkeys, "geom")

            # 2. => {TS zkeys: reagent gkeys} (map TS gkeys into reagent gkeys)
            zc_ = zmat_conv.relabel(zc_, rgt_gkey_dct, "geom")

            # 3. => {TS zkeys: reagent zkeys} (get TS => reagent zkey mapping)
            tz_to_rz = zmat_conv.subset_isomorphism(zc_, rgt_zc)

            # 4. => [reagent TS zkeys] (get sorted TS zkeys for the reagent)
            rgt_zkeys = dict_.keys_sorted_by_value(tz_to_rz)

            # 5. Save the result
            rgts_zkeys.append(rgt_zkeys)

        return rgts_zkeys

    rcts_zkeys = _get_reagent_zmat_keys(
        rgts_gkeys=reactants_keys(rxn),
        rgt_zcs=rct_zcs,
        rgt_gkey_dcts=reactant_mappings(rxn, rev=True),
    )

    prds_zkeys = _get_reagent_zmat_keys(
        rgts_gkeys=products_keys(rxn),
        rgt_zcs=prd_zcs,
        rgt_gkey_dcts=product_mappings(rxn, rev=True),
    )

    return from_data(
        tsg=ztsg,
        rcts_keys=rcts_zkeys,
        prds_keys=prds_zkeys,
        cla=class_(rxn),
        struc_typ="zmat",
        ts_zc=ts_zc,
        rct_zcs=rct_zcs,
        prd_zcs=prd_zcs,
    )


def undo_zmatrix_conversion(
    rxn: Reaction,
    ts_zc: ZmatConv = None,
    rct_zcs: list[ZmatConv] = None,
    prd_zcs: list[ZmatConv] = None,
    keep_info: bool = False,
) -> Reaction:
    """Undo a z-matrix conversion (dummy insertion + reordering) to the reaction

    This can be used to match the original geometry after conversio to z-matrix

    :param rxn: A converted reaction object
    :type rxn: Reaction
    :param ts_zc: A z-matrix conversion for the TS graph, defaults to None
    :type ts_zc: ZmatConv
    :param rct_zcs: Z-matrix conversions for reactant structures, defaults to None
    :type rct_zcs: list[ZmatConv], optional
    :param prd_zcs: Z-matrix conversions for product structures, defaults to None
    :type prd_zcs: list[ZmatConv], optional
    :param keep_zc: Keep structural conversion info?, defaults to True
    :type keep_zc: bool, optional
    :returns: The original reaction object
    :rtype: Reaction
    """
    assert not has_structures(
        rxn, complete=False
    ), f"This should not be done for a reaction with structures:\n{rxn}"

    # Read in pre-existing z-matrix conversion info, if there is some
    ts_zc = ts_conversion_info(rxn) if ts_zc is None else ts_zc
    rct_zcs = reactants_conversion_info(rxn) if rct_zcs is None else rct_zcs
    prd_zcs = products_conversion_info(rxn) if prd_zcs is None else prd_zcs

    # If the conversion info is still missing, get it from the graphs
    tsg = ts_graph(rxn)
    rct_gras = reactant_graphs(rxn, shift_keys=True)
    prd_gras = product_graphs(rxn, shift_keys=True)

    zmat_conv_ = graph.zmatrix_conversion_info
    ts_zc = zmat_conv_(tsg) if ts_zc is None else ts_zc
    rct_zcs = list(map(zmat_conv_, rct_gras)) if rct_zcs is None else rct_zcs
    prd_zcs = list(map(zmat_conv_, prd_gras)) if prd_zcs is None else prd_zcs

    # 0. Read in the current TS graph
    ztsg = ts_graph(rxn)

    # 1. Transform the TS graph
    gtsg = graph.undo_zmatrix_conversion(ztsg, ts_zc)

    # 2. Transform the reagent keys
    def _get_reagent_geom_keys(rgts_zkeys, rgt_zcs, rgt_zkey_dcts):
        """
        :param rgts_zkeys: Keys for each reagent realtive to the TS z-matrix
        :param rgt_zcs: Z-matrix conversions for each reagent
        :param rgt_zkey_dcts: Mappings from TS zkeys onto each reagent's zkeys
        """
        rgts_gkeys = []
        for rgt_zkeys, rgt_zc, rgt_zkey_dct in zip(rgts_zkeys, rgt_zcs, rgt_zkey_dcts):
            # 1. {TS zkeys: TS gkeys} (extract the subset matching reagent zkeys)
            zc_ = zmat_conv.subset(ts_zc, rgt_zkeys, "zmat")

            # 2. => {reagent zkeys: TS gkeys} (map TS zkeys into reagent zkeys)
            zc_ = zmat_conv.relabel(zc_, rgt_zkey_dct, "zmat")

            # 3. {TS gkeys: reagent gkeys} (get TS => reagent gkeys mapping)
            tg_to_rz = zmat_conv.relabel_dict(zc_, "geom")
            rz_to_rg = zmat_conv.relabel_dict(rgt_zc, "zmat")
            tg_to_rg = dict_.compose(rz_to_rg, tg_to_rz)

            # 4. [reagent TS gkeys] (get sorted TS gkeys for the reagent)
            rgt_gkeys = dict_.keys_sorted_by_value(tg_to_rg)

            # 5. Save the result
            rgts_gkeys.append(rgt_gkeys)

        return rgts_gkeys

    rcts_gkeys = _get_reagent_geom_keys(
        rgts_zkeys=reactants_keys(rxn),
        rgt_zkey_dcts=reactant_mappings(rxn, rev=True),
        rgt_zcs=rct_zcs,
    )

    prds_gkeys = _get_reagent_geom_keys(
        rgts_zkeys=products_keys(rxn),
        rgt_zkey_dcts=product_mappings(rxn, rev=True),
        rgt_zcs=prd_zcs,
    )

    return from_data(
        tsg=gtsg,
        rcts_keys=rcts_gkeys,
        prds_keys=prds_gkeys,
        cla=class_(rxn),
        struc_typ="geom" if keep_info else None,
        ts_zc=ts_zc if keep_info else None,
        rct_zcs=rct_zcs if keep_info else None,
        prd_zcs=prd_zcs if keep_info else None,
    )


def without_dummy_atoms(rxn: Reaction) -> Reaction:
    """remove dummy atoms from the reactants or products"""
    return undo_zmatrix_conversion(rxn, keep_info=False)


def remove_atoms(rxn: Reaction, keys: list[int], stereo: bool = False) -> Reaction:
    """Remove atoms from the reaction object.

    :param rxn: A reaction object
    :type rxn: Reaction
    :param keys: The keys to remove
    :type keys: Sequence[int]
    :return: The reaction object, with the atoms removed
    :rtype: Reaction
    """
    # First, relabel, moving keys to end
    rxn0 = copy.copy(rxn)
    tsg0 = ts_graph(rxn0)
    all_keys = [*sorted(graph.atom_keys(tsg0) - set(keys)), *keys]
    key_dct = dict(map(reversed, enumerate(all_keys)))
    keys = list(map(key_dct.get, keys))
    rxn = relabel(rxn, key_dct, struc=True)

    tsg = graph.remove_atoms(ts_graph(rxn), keys, stereo=stereo)
    rcts_keys = [[k for k in ks if k not in keys] for ks in reactants_keys(rxn)]
    prds_keys = [[k for k in ks if k not in keys] for ks in products_keys(rxn)]

    if structure_type(rxn) == "zmat":
        raise NotImplementedError("Cannot remove atoms from z-matrix structures")

    ts_struc = ts_structure(rxn)
    if ts_struc is not None:
        ts_struc = geom.remove(ts_struc, keys)

    rct_strucs = reactant_structures(rxn)
    if rct_strucs is not None:
        key_dcts = reactant_mappings(rxn, rev=True)
        keys_lst = [[key_dct[k] for k in keys if k in key_dct] for key_dct in key_dcts]
        rct_strucs = [
            geom.remove(g, ks) for g, ks in zip(rct_strucs, keys_lst, strict=True)
        ]

    prd_strucs = product_structures(rxn)
    if prd_strucs is not None:
        key_dcts = product_mappings(rxn, rev=True)
        keys_lst = [[key_dct[k] for k in keys if k in key_dct] for key_dct in key_dcts]
        prd_strucs = [
            geom.remove(g, ks) for g, ks in zip(prd_strucs, keys_lst, strict=True)
        ]

    return from_data(
        tsg=tsg,
        rcts_keys=rcts_keys,
        prds_keys=prds_keys,
        cla=class_(rxn),
        ts_struc=ts_struc,
        rct_strucs=rct_strucs,
        prd_strucs=prd_strucs,
        struc_typ=structure_type(rxn),
    )


def unique(rxns: list[Reaction]) -> list[Reaction]:
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


def filter_viable_reactions(rxns: list[Reaction]) -> list[Reaction]:
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

    def _sing_spin_products(rxn):
        prd_gras = product_graphs(rxn)
        mult = sum(
            map(
                graph.maximum_spin_multiplicity,
                map(graph.kekule, prd_gras),
            )
        )
        return mult == 1

    for rxn in all_rxns:
        # Check for separated radical sites
        sep_rad = _produces_separated_radical_sites(rxn)
        hi_spin = _high_spin_products(rxn)
        sing_spin = _sing_spin_products(rxn)

        # Add more conditions here, as needed ...

        if not ((sep_rad and sing_spin) or hi_spin):
            rxns.append(rxn)

    return tuple(rxns)


# helpers
def _identify_sequence_structure_type(strucs):
    """Identify the type of a sequence of structures

    :param strucs: The structures (geom or zmat)
    :type strucs: list[automol geom or zmat data structure]
    """
    strucs = [None] if strucs is None else strucs
    return (
        "geom"
        if all(map(geom.is_valid, strucs))
        else "zmat" if all(map(zmat.is_valid, strucs)) else None
    )
