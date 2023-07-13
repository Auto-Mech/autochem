""" reaction classifiers and reaction-class-specific functions

Function arguments:
    Each function takes a list of reactant graphs and a list of product graphs.
    Note that the reactant graphs *cannot* have overlapping atom keys, and
    likewise for the product graphs. Otherwise, there would be no way to
    express the bonds broken and formed between reactants.
"""

import itertools
import yaml
import numpy
import automol.geom.ts
import automol.graph
from automol import par
from automol.graph import ts


class Reaction:
    """ Describes a specific reaction

    :param class_: the name of the reaction class
    :type class_: str
    :param forward_ts_graph: a graph representing the transition state in the
        forward direction; keys must match `reactants_keys`
    :param backward_ts_graph: a graph representing the transition state in the
        backward direction; keys must match `products_keys`
    :param reactants_keys: a sequence of keys, one for each reactant, for
        extracting the reactants in order from the `forward_ts_graph`
    :type reactants_keys: tuple[tuple[int]]
    :param products_keys: a sequence of keys, one for each product, for
        extracting the products in order from the `product_ts_graph`
    :type products_keys: tuple[tuple[int]]
    """

    def __init__(self, rxn_cls, forw_tsg, back_tsg, rcts_keys, prds_keys):
        """ constructor
        """
        rcts_keys = tuple(map(tuple, map(sorted, rcts_keys)))
        prds_keys = tuple(map(tuple, map(sorted, prds_keys)))

        # Check the reaction class
        assert par.is_reaction_class(rxn_cls), (
            f"{rxn_cls} is not a reaction class")

        # Check the reactant keys and the forward transition state graph
        all_rcts_keys = set(itertools.chain(*rcts_keys))
        forw_keys = automol.graph.atom_keys(forw_tsg)
        assert all_rcts_keys == forw_keys, (
            f"{str(all_rcts_keys)} != {str(forw_keys)}")

        # Check the product keys and the backward transition state graph
        all_prds_keys = set(itertools.chain(*prds_keys))
        back_keys = automol.graph.atom_keys(back_tsg)
        assert all_prds_keys == back_keys, (
            f"{str(all_prds_keys)} != {str(back_keys)}")

        # Check that the stereo information is consistent
        # (All or nothing - no partially assigned stereo)
        if automol.graph.has_stereo(forw_tsg):
            ste_keys = automol.graph.stereogenic_keys(forw_tsg)
            assert not ste_keys, (
                f"TS graph has unassigned stereo at {ste_keys}:\n{forw_tsg}")

        # Check that the forward and reverse TS graphs are consistent
        forw_tsg_comp = automol.graph.without_dummy_atoms(forw_tsg)
        back_tsg_comp = automol.graph.without_dummy_atoms(back_tsg)
        assert automol.graph.isomorphism(
            ts.reverse(forw_tsg_comp), back_tsg_comp)

        # Set attributes
        self.class_ = rxn_cls
        self.reactants_keys = rcts_keys
        self.products_keys = prds_keys
        self.forward_ts_graph = forw_tsg
        self.backward_ts_graph = back_tsg

    def copy(self):
        """ return a copy of this Reaction
        """
        return Reaction(
            self.class_, self.forward_ts_graph, self.backward_ts_graph,
            self.reactants_keys, self.products_keys)

    def __eq__(self, other):
        """ equality operator
        """
        if not isinstance(other, Reaction):
            ret = False
        else:
            ret = (self.class_ == other.class_ and
                   self.forward_ts_graph == other.forward_ts_graph and
                   self.backward_ts_graph == other.backward_ts_graph and
                   self.reactants_keys == other.reactants_keys and
                   self.products_keys == other.products_keys)
        return ret

    def __repr__(self):
        """ string representation of the object
        """
        return string(self)


def string(rxn, one_indexed=True):
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

    forw_ts_dct = automol.graph.yaml_dictionary(rxn.forward_ts_graph,
                                                one_indexed=one_indexed)
    back_ts_dct = automol.graph.yaml_dictionary(rxn.backward_ts_graph,
                                                one_indexed=one_indexed)
    yaml_dct = {}
    yaml_dct['reaction class'] = rxn.class_
    yaml_dct['forward TS atoms'] = forw_ts_dct['atoms']
    yaml_dct['forward TS bonds'] = forw_ts_dct['bonds']
    yaml_dct['reactants keys'] = rcts_keys
    yaml_dct['backward TS atoms'] = back_ts_dct['atoms']
    yaml_dct['backward TS bonds'] = back_ts_dct['bonds']
    yaml_dct['products keys'] = prds_keys
    rxn_str = yaml.dump(
        yaml_dct, default_flow_style=None, sort_keys=False)
    return rxn_str


def from_string(rxn_str, one_indexed=True):
    """ Write a reaction object to a string

    :param rxn_str: string containing the reaction object
    :type rxn_str: str
    :param one_indexed: parameter to store keys in one-indexing
    :type one_indexed: bool
    :rtype: Reaction
    """
    yaml_dct = yaml.load(rxn_str, Loader=yaml.FullLoader)

    rxn_cls = yaml_dct['reaction class']
    rcts_keys = yaml_dct['reactants keys']
    prds_keys = yaml_dct['products keys']

    if one_indexed:
        rcts_keys = [[k-1 for k in ks] for ks in rcts_keys]
        prds_keys = [[k-1 for k in ks] for ks in prds_keys]

    forw_ts_dct = {
        'atoms': yaml_dct['forward TS atoms'],
        'bonds': yaml_dct['forward TS bonds']}
    forw_tsg = automol.graph.from_yaml_dictionary(forw_ts_dct,
                                                  one_indexed=one_indexed)
    back_ts_dct = {
        'atoms': yaml_dct['backward TS atoms'],
        'bonds': yaml_dct['backward TS bonds']}
    back_tsg = automol.graph.from_yaml_dictionary(back_ts_dct,
                                                  one_indexed=one_indexed)
    rxn = Reaction(rxn_cls, forw_tsg, back_tsg, rcts_keys, prds_keys)
    return rxn


def reverse(rxn):
    """ Obtains the reaction object for the reverse reaction

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: Reaction
    """
    rxn_cls = par.reverse_reaction_class(rxn.class_)
    if rxn_cls is not None:
        forw_tsg = rxn.backward_ts_graph
        back_tsg = rxn.forward_ts_graph
        rcts_keys = rxn.products_keys
        prds_keys = rxn.reactants_keys
        rxn = Reaction(rxn_cls, forw_tsg, back_tsg, rcts_keys, prds_keys)
    else:
        rxn = None
    return rxn


def atom_mapping(rxn, rev=False, stereo=True):
    """ Determine a mapping from reaction atoms to product atoms.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle reaction direction
    :type rev: bool
    :param stereo: Check for stereo parities in the mapping?
    :type stereo: bool
    :returns: the mapping from reactant atoms to product atoms
    :rtype: dict
    """
    iso_dct = automol.graph.isomorphism(
        ts.reverse(rxn.forward_ts_graph), rxn.backward_ts_graph,
        stereo=stereo, dummy=False)
    if rev:
        iso_dct = dict(map(reversed, iso_dct.items()))

    return iso_dct


def sort_order(rxn):
    """ Get the sort order for reactants and products, based on their keys

    :param rxn: the reaction object
    :type rxn: Reaction
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


def forming_bond_keys(rxn, rev=False):
    """ Obtain forming bonds for the reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle reaction direction
    :type rev: bool
    :rtype: frozenset[frozenset[int]]
    """
    if rev:
        tsg = rxn.backward_ts_graph
    else:
        tsg = rxn.forward_ts_graph
    return ts.ts_forming_bond_keys(tsg)


def breaking_bond_keys(rxn, rev=False):
    """ Obtain breaking bonds for the reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle reaction direction
    :type rev: bool
    :rtype: frozenset[frozenset[int]]
    """
    if rev:
        tsg = rxn.backward_ts_graph
    else:
        tsg = rxn.forward_ts_graph
    return ts.ts_breaking_bond_keys(tsg)


def forming_rings_atom_keys(rxn, rev=False):
    """ Obtain atom keys for rings containing at least one forming bond.

    Atom keys are sorted by connectivity.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle reaction direction
    :type rev: bool
    :rtype: tuple[tuple[int]]
    """
    if rev:
        tsg = rxn.backward_ts_graph
    else:
        tsg = rxn.forward_ts_graph
    return ts.forming_rings_atom_keys(tsg)


def forming_rings_bond_keys(rxn, rev=False):
    """ Obtain bond keys for rings containing at least one forming bond.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle reaction direction
    :type rev: bool
    :rtype: tuple[frozenset[frozenset[int]]]
    """
    if rev:
        tsg = rxn.backward_ts_graph
    else:
        tsg = rxn.forward_ts_graph
    return ts.forming_rings_bond_keys(tsg)


def breaking_rings_atom_keys(rxn, rev=False):
    """ Obtain atom keys for rings containing at least one breaking bond.

    Atom keys are sorted by connectivity.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle reaction direction
    :type rev: bool
    :rtype: tuple[tuple[int]]
    """
    if rev:
        tsg = rxn.backward_ts_graph
    else:
        tsg = rxn.forward_ts_graph
    return ts.breaking_rings_atom_keys(tsg)


def breaking_rings_bond_keys(rxn, rev=False):
    """ Obtain bond keys for rings containing at least one breaking bond.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle reaction direction
    :type rev: bool
    :rtype: tuple[frozenset[frozenset[int]]]
    """
    if rev:
        tsg = rxn.backward_ts_graph
    else:
        tsg = rxn.forward_ts_graph
    return ts.breaking_rings_bond_keys(tsg)


def reactant_graphs(rxn, rev=False):
    """ Obtain graphs of the reactants in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: tuple of automol graph data structures
    """
    if rev:
        rcts_gra = ts.products_graph(rxn.forward_ts_graph)
    else:
        rcts_gra = ts.reactants_graph(rxn.forward_ts_graph)
    rct_gras = [automol.graph.subgraph(rcts_gra, keys, stereo=True)
                for keys in rxn.reactants_keys]
    return tuple(rct_gras)


def product_graphs(rxn):
    """ Obtain graphs of the products in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: tuple of automol graph data structures
    """
    prds_gra = ts.reactants_graph(rxn.backward_ts_graph)
    prd_gras = [automol.graph.subgraph(prds_gra, keys, stereo=True)
                for keys in rxn.products_keys]
    return tuple(prd_gras)


def reactants_graph(rxn, rev=False):
    """ Obtain a (single) graph of the reactants in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle reaction direction
    :type rev: bool
    :rtype: automol graph data structure
    """
    if rev:
        tsg = rxn.backward_ts_graph
    else:
        tsg = rxn.forward_ts_graph
    return ts.reactants_graph(tsg)


def products_graph(rxn):
    """ Obtain a (single) graph of the products in this reaction.

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: automol graph data structure
    """
    tsg = rxn.backward_ts_graph
    return ts.reactants_graph(tsg)


def standard_keys(rxn):
    """ Standardize keys for the reaction object

    Ensures that keys follow zero-indexing in order with no skips

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: A copy of the reaction object, with standardized keys
    :rtype: Reaction
    """
    rct_keys = list(map(sorted, rxn.reactants_keys))
    prd_keys = list(map(sorted, rxn.products_keys))
    rct_key_dct = {k: i for i, k in enumerate(itertools.chain(*rct_keys))}
    prd_key_dct = {k: i for i, k in enumerate(itertools.chain(*prd_keys))}

    rxn_cls = rxn.class_
    forw_tsg = automol.graph.relabel(rxn.forward_ts_graph, rct_key_dct)
    back_tsg = automol.graph.relabel(rxn.backward_ts_graph, prd_key_dct)
    rcts_keys = tuple(tuple(map(rct_key_dct.__getitem__, keys))
                      for keys in rxn.reactants_keys)
    prds_keys = tuple(tuple(map(prd_key_dct.__getitem__, keys))
                      for keys in rxn.products_keys)
    rxn = Reaction(rxn_cls, forw_tsg, back_tsg, rcts_keys, prds_keys)
    return rxn


def has_standard_keys(rxn):
    """ Determine if this reaction has standard keys

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return rxn == standard_keys(rxn)


def standard_keys_with_sorted_geometries(rxn, rct_geos, prd_geos):
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


def relabel(rxn, key_dct, prod=False):
    """ Relabel keys in the reactants or products

    :param rxn: the reaction object
    :type rxn: Reaction
    :param key_dct: A dictionary mapping current keys into new keys
    :type key_dct: dict[int: int]
    :param prod: Apply this operation to the products instead of the reactants
    :type prod: bool
    :returns: A relabeled reaction object
    :rtype: Reaction
    """
    if prod:
        tsg = rxn.backward_ts_graph
        keys_lst = rxn.products_keys
    else:
        tsg = rxn.forward_ts_graph
        keys_lst = rxn.reactants_keys

    tsg = automol.graph.relabel(tsg, key_dct)
    keys_lst = [sorted(map(key_dct.__getitem__, keys)) for keys in keys_lst]

    rxn_cls = rxn.class_
    forw_tsg = rxn.forward_ts_graph
    back_tsg = rxn.backward_ts_graph
    rcts_keys = rxn.reactants_keys
    prds_keys = rxn.products_keys
    if prod:
        back_tsg = tsg
        prds_keys = keys_lst
    else:
        forw_tsg = tsg
        rcts_keys = keys_lst
    rxn = Reaction(rxn_cls, forw_tsg, back_tsg, rcts_keys, prds_keys)
    return rxn


def without_stereo(rxn):
    """ Remove all stereo information from the reaction object

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the reaction object, without stereo information
    :rtype: Reaction
    """
    rxn_cls = rxn.class_
    forw_tsg = automol.graph.without_stereo(rxn.forward_ts_graph)
    back_tsg = automol.graph.without_stereo(rxn.backward_ts_graph)
    rcts_keys = rxn.reactants_keys
    prds_keys = rxn.products_keys
    rxn = Reaction(rxn_cls, forw_tsg, back_tsg, rcts_keys, prds_keys)
    return rxn


def add_dummy_atoms(rxn, dummy_key_dct, product=False):
    """ Add dummy atoms to the reactants or products

    :param rxn: the reaction object
    :type rxn: Reaction
    :param dummy_key_dct: Keys of new dummy atoms, by atom that they connect to
    :type dummy_key_dct: dict[int: int]
    :param prod: Apply this operation to the products instead of the reactants
    :type prod: bool
    :returns: A reaction object with dummy atoms added
    :rtype: Reaction
    """
    if product:
        tsg = rxn.backward_ts_graph
        keys_lst = rxn.products_keys
    else:
        tsg = rxn.forward_ts_graph
        keys_lst = rxn.reactants_keys

    tsg = automol.graph.add_dummy_atoms(tsg, dummy_key_dct)
    keys_lst = list(map(list, keys_lst))
    for key, dummy_key in dummy_key_dct.items():
        for keys in keys_lst:
            if key in keys:
                keys.append(dummy_key)
                break

    rxn_cls = rxn.class_
    forw_tsg = rxn.forward_ts_graph
    back_tsg = rxn.backward_ts_graph
    rcts_keys = rxn.reactants_keys
    prds_keys = rxn.products_keys
    if product:
        back_tsg = tsg
        prds_keys = keys_lst
    else:
        forw_tsg = tsg
        rcts_keys = keys_lst
    rxn = Reaction(rxn_cls, forw_tsg, back_tsg, rcts_keys, prds_keys)
    return rxn


def insert_dummy_atoms(rxn, dummy_key_dct, product=False):
    """ insert dummy atoms into the reactants or products

    :param dummy_key_dct: dummy atom keys, by key of the atom they are
        connected to
    :param product: insert dummy atoms into the products instead of the
        reactants?
    """
    for key, dummy_key in reversed(sorted(dummy_key_dct.items())):
        rxn = _insert_dummy_atom(rxn, key, dummy_key, product=product)

    return rxn


def _insert_dummy_atom(rxn, key, dummy_key, product=False):
    tsg = rxn.backward_ts_graph if product else rxn.forward_ts_graph
    keys = sorted(automol.graph.atom_keys(tsg))
    dummy_key_ = max(keys) + 1
    rxn = add_dummy_atoms(rxn, {key: dummy_key_}, product=product)

    if dummy_key != dummy_key_:
        assert dummy_key in keys
        idx = keys.index(dummy_key)
        key_dct = {}
        key_dct.update({k: k for k in keys[:idx]})
        key_dct[dummy_key_] = dummy_key
        key_dct.update({k: k+1 for k in keys[idx:]})

    rxn = relabel(rxn, key_dct)
    return rxn


def without_dummy_atoms(rxn, product=False):
    """ remove dummy atoms from the reactants or products

    :param product: do this do the products instead of the reactants?
    """
    if product:
        tsg = rxn.backward_ts_graph
        keys_lst = rxn.products_keys
    else:
        tsg = rxn.forward_ts_graph
        keys_lst = rxn.reactants_keys

    dummy_keys = automol.graph.atom_keys(tsg, symb='X')

    tsg = automol.graph.remove_atoms(tsg, dummy_keys, stereo=True)
    keys_lst = [[k for k in ks if k not in dummy_keys] for ks in keys_lst]

    rxn_cls = rxn.class_
    forw_tsg = rxn.forward_ts_graph
    back_tsg = rxn.backward_ts_graph
    rcts_keys = rxn.reactants_keys
    prds_keys = rxn.products_keys
    if product:
        back_tsg = tsg
        prds_keys = keys_lst
    else:
        forw_tsg = tsg
        rcts_keys = keys_lst
    rxn = Reaction(rxn_cls, forw_tsg, back_tsg, rcts_keys, prds_keys)
    return rxn


def relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct, product=False):
    """ relabel the reaction object to correspond with a z-matrix converted
    from a geometry

    :param rxn: the reaction object
    :param zma_keys: graph keys in the order they appear in the z-matrix
    :param dummy_key_dct: dummy keys introduced on z-matrix conversion, by atom
        they are attached to
    :param product: do this for the products instead of the reactants?
    :type product: bool
    """
    rxn = add_dummy_atoms(rxn, dummy_key_dct, product=product)
    key_dct = dict(map(reversed, enumerate(zma_keys)))
    rxn = relabel(rxn, key_dct, prod=product)
    return rxn


def relabel_for_geometry(rxn, product=False):
    """ relabel the reaction object to correspond with a geometry converted
    from a z-matrix

    :param rxn: the reaction object
    :param product: do this do the products instead of the reactants?
    """
    rxn = without_dummy_atoms(rxn, product=product)

    if product:
        tsg = rxn.backward_ts_graph
    else:
        tsg = rxn.forward_ts_graph

    keys = sorted(automol.graph.atom_keys(tsg))
    key_dct = dict(map(reversed, enumerate(keys)))
    rxn = relabel(rxn, key_dct, prod=product)
    return rxn


def reaction_class(rxn):
    """return the reaction class string

    :param rxn: the reaction object
    """
    return rxn.class_


def unique(rxns):
    """ Get reactions with distinct TSs from a list with redundancies

    :param rxns: a sequence of reaction objects
    :returns: unique reaction objects
    """
    all_rxns = rxns
    rxns = []

    def isomorphic_(rxn1, rxn2):
        tsg1 = rxn1.forward_ts_graph
        tsg2 = rxn2.forward_ts_graph
        return automol.graph.isomorphic(tsg1, tsg2, stereo=True)

    for rxn in all_rxns:
        if not any(isomorphic_(rxn, r) for r in rxns):
            rxns.append(rxn)

    return tuple(rxns)


def is_radical_radical(zrxn, rev=False):
    """ Is this a radical-radical reaction

    :param rxn: the reaction object
    :type rxn: Reaction
    :param rev: parameter to toggle reaction direction
    :type rev: bool
    :rtype: boolean
    """
    _is_rad_rad = False
    tsg = reactant_graphs(zrxn, rev=rev)
    if len(tsg) == 2:
        rct_i, rct_j = tsg
        if (automol.graph.is_radical_species(rct_i)
                and automol.graph.is_radical_species(rct_j)):
            _is_rad_rad = True
    return _is_rad_rad


def is_barrierless(zrxn, rev=False):
    """ Is this a barrierless reaction

        :param rxn: the reaction object
        :type rxn: Reaction
        :param rev: parameter to toggle reaction direction
        :type rev: bool
        :rtype: boolean
    """
    _is_barr = False
    if is_radical_radical(zrxn, rev=rev):
        _is_barr = True
    # Add other examples
    return _is_barr


def filter_viable_reactions(rxns):
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
