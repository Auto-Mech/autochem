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
from typing import List
import yaml
import numpy
import automol.geom.ts
import automol.graph
from automol import par
from automol.graph import ts


@dataclasses.dataclass
class Reaction:
    """ Describes a specific reaction

    :param class_: The reaction class
    :type class_: str
    :param ts_graph: The transition state graph
    :type ts_graph: automol graph data structure
    :param reactants_keys: Keys to the reactant molecules in `tsg`
    :type reactants_keys: tuple[tuple[int]]
    :param products_keys: Keys to the product molecules in `tsg`
    :type products_keys: tuple[tuple[int]]
    """
    class_: str
    ts_graph: tuple
    reactants_keys: tuple
    products_keys: tuple

    def __repr__(self):
        """ string representation of the object
        """
        return string(self)


def from_forward_reverse(cla, ftsg, rtsg, rcts_keys, prds_keys):
    """ Construct a Reaction dataclass from forward and reverse TS graphs

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
    """
    # Determine the reaction mapping
    rmap_dct = automol.graph.isomorphism(ts.reverse(rtsg), ftsg, dummy=False)
    rcts_keys = list(map(sorted, rcts_keys))
    prds_keys = list(map(sorted, prds_keys))

    # Reverse-map the products keys so they line up with the forward TS graph
    prds_keys = tuple(tuple(map(rmap_dct.__getitem__, ks)) for ks in prds_keys)

    return from_data(cla, ftsg, rcts_keys, prds_keys)


def from_data(cla, tsg, rcts_keys, prds_keys):
    """ Construct a Reaction dataclass from data

    :param cla: The reaction class
    :type cla: str
    :param tsg: The transition state graph
    :type tsg: automol graph data structure
    :param rcts_keys: Keys to the reactant molecules in `tsg`
    :type rcts_keys: tuple[tuple[int]]
    :param prds_keys: Keys to the product molecules in `tsg`
    :type prds_keys: tuple[tuple[int]]
    :param rmap: A dictionary mapping reactant atoms onto product atoms
    :type rmap: dict[int: int]
    """

    # Check the reaction class
    assert par.is_reaction_class(cla), (
        f"{cla} is not a reaction class")

    # Check the TS graph...
    # If present, stereo information should be complete
    if automol.graph.has_stereo(tsg):
        ste_keys = automol.graph.stereogenic_keys(tsg)
        assert not ste_keys, (
            f"TS graph has unassigned stereo at {ste_keys}:\n{tsg}")

    # Check the reactants and products keys
    rcts_keys = tuple(map(tuple, rcts_keys))
    prds_keys = tuple(map(tuple, prds_keys))

    return Reaction(
        class_=cla, ts_graph=tsg, reactants_keys=rcts_keys,
        products_keys=prds_keys)


def from_string(rxn_str, one_indexed=True):
    """ Write a reaction object to a string

    :param rxn_str: string containing the reaction object
    :type rxn_str: str
    :param one_indexed: parameter to store keys in one-indexing
    :type one_indexed: bool
    :rtype: Reaction
    """
    yaml_dct = yaml.load(rxn_str, Loader=yaml.FullLoader)

    cla = yaml_dct['reaction class']
    tsg = automol.graph.from_yaml_dictionary(yaml_dct, one_indexed=one_indexed)
    rcts_keys = tuple(map(tuple, yaml_dct['reactants keys']))
    prds_keys = tuple(map(tuple, yaml_dct['products keys']))

    if one_indexed:
        rcts_keys = [[k-1 for k in ks] for ks in rcts_keys]
        prds_keys = [[k-1 for k in ks] for ks in prds_keys]

    return from_data(cla=cla, tsg=tsg, rcts_keys=rcts_keys,
                     prds_keys=prds_keys)


def string(rxn: Reaction, one_indexed=True):
    """ Write a reaction object to a string

    :param rxn: the reaction object
    :type rxn: Reaction
    :param one_indexed: parameter to store keys in one-indexing
    :type one_indexed: bool
    :rtype: str
    """
    rcts_keys = list(map(list, rxn.reactants_keys))
    prds_keys = list(map(list, rxn.products_keys))

    if one_indexed:
        rcts_keys = [[k+1 for k in ks] for ks in rcts_keys]
        prds_keys = [[k+1 for k in ks] for ks in prds_keys]

    yaml_dct = automol.graph.yaml_dictionary(
        rxn.ts_graph, one_indexed=one_indexed)
    yaml_dct['reactants keys'] = list(map(list, rcts_keys))
    yaml_dct['products keys'] = list(map(list, prds_keys))
    yaml_dct['reaction class'] = rxn.class_
    rxn_str = yaml.dump(
        yaml_dct, default_flow_style=None, sort_keys=False)
    return rxn_str


def relabel(rxn: Reaction, key_dct):
    """ Relabel keys in the reactants or products

    :param rxn: the reaction object
    :type rxn: Reaction
    :param key_dct: A dictionary mapping current keys into new keys
    :type key_dct: dict[int: int]
    :returns: A relabeled reaction object
    :rtype: Reaction
    """
    rxn = copy.deepcopy(rxn)
    rxn.ts_graph = automol.graph.relabel(rxn.ts_graph, key_dct)
    rxn.reactants_keys = tuple(tuple(map(key_dct.__getitem__, ks))
                               for ks in rxn.reactants_keys)
    rxn.products_keys = tuple(tuple(map(key_dct.__getitem__, ks))
                              for ks in rxn.products_keys)
    return rxn


def reverse(rxn: Reaction):
    """ Obtains the reaction object for the reverse reaction

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: Reaction
    """
    return from_data(
        cla=par.reverse_reaction_class(rxn.class_),
        tsg=automol.graph.ts.reverse(rxn.ts_graph),
        rcts_keys=rxn.products_keys,
        prds_keys=rxn.reactants_keys,
    )


def mapping(rxn: Reaction, inp, out):
    """ Determine a mapping from TS atoms to reactant atoms.

    Code:
        'T': TS keys
        'R': reactant keys
        'P': product keys

    :param rxn: the reaction object
    :type rxn: Reaction
    :param inp: which keys will go into the mapping, 'T', 'R', or 'P'
    :type inp: str
    :param out: which keys will come out of the mapping, 'T', 'R', or 'P'
    :type out: str
    :returns: A dictionary mapping `inp` keys into `out` keys
    :rtype: dict
    """
    keys = sorted(automol.graph.atom_keys(rxn.ts_graph))
    rct_keys = list(itertools.chain(*rxn.reactants_keys))
    prd_keys = list(itertools.chain(*rxn.products_keys))

    key_dct = {
        'T': keys,
        'R': [rct_keys.index(k) for k in keys],
        'P': [prd_keys.index(k) for k in keys],
    }
    inp_keys = key_dct[inp]
    out_keys = key_dct[out]
    map_dct = dict(zip(inp_keys, out_keys))
    return map_dct


def reaction_mapping(rxn: Reaction, rev=False):
    """ Determine a mapping from reactant atoms to product atoms.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle direction of the mapping
    :type rev: bool
    :returns: The mapping
    :rtype: dict
    """
    rct_keys = list(itertools.chain(*rxn.reactants_keys))
    prd_keys = list(itertools.chain(*rxn.products_keys))
    map_dct = (dict(zip(prd_keys, rct_keys)) if rev else
               dict(zip(rct_keys, rct_keys)))
    return map_dct


def sort_order(rxn: Reaction):
    """ Get the sort order for reactants and products, based on their keys

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle reaction direction
    :type rev: bool
    :returns: The sort order of the reactants and products
    :rtype: (tuple[int], tuple[int])
    """
    if len(rxn.reactants_keys) == 1:
        rct_idxs = [0]
    else:
        rct_keys = numpy.empty((len(rxn.reactants_keys),), dtype=object)
        rct_keys[:] = list(map(tuple, map(sorted, rxn.reactants_keys)))
        rct_idxs = numpy.argsort(rct_keys)

    if len(rxn.products_keys) == 1:
        prd_idxs = [0]
    else:
        prd_keys = numpy.empty((len(rxn.products_keys),), dtype=object)
        prd_keys[:] = list(map(tuple, map(sorted, rxn.products_keys)))
        prd_idxs = numpy.argsort(prd_keys)

    rct_idxs, prd_idxs = map(tuple, (rct_idxs, prd_idxs))
    return rct_idxs, prd_idxs


def reactant_graphs(rxn: Reaction, shift_keys=False):
    """ Obtain graphs of the reactants in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param shift_keys: Shift keys after first reagent, to prevent overlap?
    :type shift_keys: bool
    :param overlap: Return zero-indexed 
    :rtype: tuple of automol graph data structures
    """
    map_dct = mapping(rxn, 'T', 'R')

    rcts_gra = ts.reactants_graph(rxn.ts_graph)
    rct_gras = [automol.graph.subgraph(rcts_gra, ks, stereo=True)
                for ks in rxn.reactants_keys]
    rct_gras = [automol.graph.relabel(g, map_dct, check=False)
                for g in rct_gras]
    if not shift_keys:
        rct_gras = [automol.graph.standard_keys(g) for g in rct_gras]
    return tuple(rct_gras)


def product_graphs(rxn: Reaction, shift_keys=False):
    """ Obtain graphs of the products in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param shift_keys: Shift keys after first reagent, to prevent overlap?
    :type shift_keys: bool
    :rtype: tuple of automol graph data structures
    """
    return reactant_graphs(reverse(rxn), shift_keys=shift_keys)


def reactants_graph(rxn: Reaction):
    """ Obtain a (single) graph of the reactants in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: automol graph data structure
    """
    map_dct = mapping(rxn, 'T', 'R')
    rcts_gra = ts.reactants_graph(rxn.ts_graph)
    return automol.graph.relabel(rcts_gra, map_dct, check=True)


def products_graph(rxn: Reaction):
    """ Obtain a (single) graph of the products in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: automol graph data structure
    """
    return reactants_graph(reverse(rxn))


def standard_keys(rxn: Reaction):
    """ Standardize keys for the reaction object

    Ensures that keys follow zero-indexing in order with no skips

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: A copy of the reaction object, with standardized keys
    :rtype: Reaction
    """
    rct_keys = list(map(sorted, rxn.reactants_keys))
    rct_key_dct = {k: i for i, k in enumerate(itertools.chain(*rct_keys))}

    rxn = relabel(rxn, rct_key_dct)
    return rxn


def has_standard_keys(rxn: Reaction):
    """ Determine if this reaction has standard keys

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return rxn == standard_keys(rxn)


def standard_keys_with_sorted_geometries(rxn: Reaction, rct_geos, prd_geos):
    """ Standardize keys and line up geometries to match proper
    reactant/product sort ordering for the reaction

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rct_geos: Geometries for the reactants, with indices matching `rxn`
    :type rct_geos: tuple of automol geometry objects
    :param prd_geos: Geometries for the products, with indices matching `rxn`
    :type prd_geos: tuple of automol geometry objects
    :returns: A standardized reaction object, with matching re-sorted reactants
        and products
    :rtype: (Reaction, automol geometry objects, automol geometry objects)
    """
    rct_idxs, prd_idxs = sort_order(rxn)
    rct_geos = tuple(map(rct_geos.__getitem__, rct_idxs))
    prd_geos = tuple(map(prd_geos.__getitem__, prd_idxs))
    rxn = standard_keys(rxn)
    return rxn, rct_geos, prd_geos


def without_stereo(rxn: Reaction):
    """ Remove all stereo information from the reaction object

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the reaction object, without stereo information
    :rtype: Reaction
    """
    rxn = copy.deepcopy(rxn)
    rxn.ts_graph = automol.graph.without_stereo(rxn.ts_graph)
    return rxn


def add_dummy_atoms(rxn: Reaction, dummy_key_dct):
    """ Add dummy atoms to the reactants or products

    :param rxn: the reaction object
    :type rxn: Reaction
    :param dummy_key_dct: Keys of new dummy atoms, by atom that they connect to
    :type dummy_key_dct: dict[int: int]
    :returns: A reaction object with dummy atoms added
    :rtype: Reaction
    """
    tsg = rxn.ts_graph
    keys_lst = rxn.reactants_keys

    rxn = copy.deepcopy(rxn)

    tsg = automol.graph.add_dummy_atoms(tsg, dummy_key_dct)
    keys_lst = list(map(list, keys_lst))
    for key, dummy_key in dummy_key_dct.items():
        for keys in keys_lst:
            if key in keys:
                keys.append(dummy_key)
                break

    rxn.ts_graph = tsg
    rxn.reactants_keys = tuple(map(tuple, keys_lst))
    return rxn


def insert_dummy_atoms(rxn: Reaction, dummy_key_dct):
    """ insert dummy atoms into the reactants or products

    :param dummy_key_dct: dummy atom keys, by key of the atom they are
        connected to
    """
    for key, dummy_key in reversed(sorted(dummy_key_dct.items())):
        rxn = _insert_dummy_atom(rxn, key, dummy_key)

    return rxn


def _insert_dummy_atom(rxn, key, dummy_key):
    tsg = rxn.ts_graph
    keys = sorted(automol.graph.atom_keys(tsg))
    dummy_key_ = max(keys) + 1
    rxn = add_dummy_atoms(rxn, {key: dummy_key_})

    if dummy_key != dummy_key_:
        assert dummy_key in keys
        idx = keys.index(dummy_key)
        key_dct = {}
        key_dct.update({k: k for k in keys[:idx]})
        key_dct[dummy_key_] = dummy_key
        key_dct.update({k: k+1 for k in keys[idx:]})

    rxn = relabel(rxn, key_dct)
    return rxn


def without_dummy_atoms(rxn: Reaction):
    """ remove dummy atoms from the reactants or products

    :param product: do this do the products instead of the reactants?
    """
    keys_lst = rxn.reactants_keys
    tsg = rxn.ts_graph

    rxn = copy.deepcopy(rxn)

    dummy_keys = automol.graph.atom_keys(tsg, symb='X')
    rxn.ts_graph = automol.graph.remove_atoms(tsg, dummy_keys, stereo=True)
    rxn.reactants_keys = [
        [k for k in ks if k not in dummy_keys] for ks in keys_lst]
    return rxn


def relabel_for_zmatrix(rxn: Reaction, zma_keys, dummy_key_dct):
    """ relabel the reaction object to correspond with a z-matrix converted
    from a geometry

    :param rxn: the reaction object
    :param zma_keys: graph keys in the order they appear in the z-matrix
    :param dummy_key_dct: dummy keys introduced on z-matrix conversion, by atom
        they are attached to
    """
    rxn = add_dummy_atoms(rxn, dummy_key_dct)
    key_dct = dict(map(reversed, enumerate(zma_keys)))
    rxn = relabel(rxn, key_dct)
    return rxn


def relabel_for_geometry(rxn: Reaction):
    """ relabel the reaction object to correspond with a geometry converted
    from a z-matrix

    :param rxn: the reaction object
    :param product: do this do the products instead of the reactants?
    """
    rxn = without_dummy_atoms(rxn)

    tsg = rxn.ts_graph
    keys = sorted(automol.graph.atom_keys(tsg))
    key_dct = dict(map(reversed, enumerate(keys)))
    rxn = relabel(rxn, key_dct)
    return rxn


def reaction_class(rxn: Reaction):
    """return the reaction class string

    :param rxn: the reaction object
    """
    return rxn.class_


def unique(rxns: List[Reaction]):
    """ Get reactions with distinct TSs from a list with redundancies

    :param rxns: a sequence of reaction objects
    :returns: unique reaction objects
    """
    all_rxns = rxns
    rxns = []

    def isomorphic_(rxn1, rxn2):
        tsg1 = rxn1.ts_graph
        tsg2 = rxn2.ts_graph
        return automol.graph.isomorphic(tsg1, tsg2, stereo=True)

    for rxn in all_rxns:
        if not any(isomorphic_(rxn, r) for r in rxns):
            rxns.append(rxn)

    return tuple(rxns)


def is_radical_radical(zrxn: Reaction):
    """ Is this a radical-radical reaction

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: boolean
    """
    _is_rad_rad = False
    rct_gras = reactant_graphs(zrxn)
    if len(rct_gras) == 2:
        rct_i, rct_j = rct_gras
        if (automol.graph.is_radical_species(rct_i)
                and automol.graph.is_radical_species(rct_j)):
            _is_rad_rad = True
    return _is_rad_rad


def filter_viable_reactions(rxns: List[Reaction]):
    """ Filter a list of reactions to only include the viable ones

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
        sep_rad = any(automol.graph.has_separated_radical_sites(prd_gra)
                      for prd_gra in prd_gras)
        return sep_rad

    def _high_spin_products(rxn):
        prd_gras = product_graphs(rxn)
        mult = sum(map(automol.graph.maximum_spin_multiplicity,
                       map(automol.graph.kekule, prd_gras)))
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
