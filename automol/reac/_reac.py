""" reaction classifiers and reaction-class-specific functions

Function arguments:
    Each function takes a list of reactant graphs and a list of product graphs.
    Note that the reactant graphs *cannot* have overlapping atom keys, and
    likewise for the product graphs. Otherwise, there would be no way to
    express the bonds broken and formed between reactants.
"""
import itertools
import numpy
import automol.convert.graph
import automol.geom.ts
from automol import par
from automol.graph import ts
from automol.graph import subgraph
from automol.graph import atom_keys
from automol.graph import full_isomorphism
from automol.graph import without_dummy_atoms
from automol.graph import add_bonded_dummy_atoms


class Reaction:
    """ Describes a specific reaction

    Methods with trailing underscores have a side-effect

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
            "{} is not a reaction class".format(rxn_cls))

        # Check the reactant keys and the forward transition state graph
        all_rcts_keys = set(itertools.chain(*rcts_keys))
        assert all_rcts_keys == atom_keys(forw_tsg), (
            "{} != {}".format(str(all_rcts_keys), str(atom_keys(forw_tsg))))

        # Check the product keys and the backward transition state graph
        all_prds_keys = set(itertools.chain(*prds_keys))
        assert all_prds_keys == atom_keys(back_tsg), (
            "{} != {}".format(str(all_prds_keys), str(atom_keys(back_tsg))))

        # Check that the reactants and products are consistent
        assert full_isomorphism(ts.reverse(forw_tsg, dummies=False),
                                without_dummy_atoms(back_tsg))

        # Set attributes
        self.class_ = rxn_cls
        self.reactants_keys = rcts_keys
        self.products_keys = prds_keys
        self.forward_ts_graph = forw_tsg
        self.backward_ts_graph = back_tsg

    def sort_order(self):
        """ determine the appropriate sort order for reactants and products,
        based on their keys
        """
        if len(self.reactants_keys) == 1:
            rct_idxs = [0]
        else:
            rct_keys = list(map(tuple, map(sorted, self.reactants_keys)))
            rct_idxs = numpy.argsort(numpy.array(rct_keys, dtype=object))

        if len(self.products_keys) == 1:
            prd_idxs = [0]
        else:
            prd_keys = list(map(tuple, map(sorted, self.products_keys)))
            prd_idxs = numpy.argsort(numpy.array(prd_keys, dtype=object))

        rct_idxs, prd_idxs = map(tuple, (rct_idxs, prd_idxs))
        return rct_idxs, prd_idxs

    def key_map(self, reverse=False):
        """ get the key map taking atoms from the reactant into atoms from the
        product
        """
        iso_dct = full_isomorphism(ts.reverse(self.forward_ts_graph),
                                   self.backward_ts_graph)
        if reverse:
            iso_dct = dict(map(reversed, iso_dct.items()))

        return iso_dct

    def reverse_(self):
        """ get the reaction class for the reverse reaction
        """
        self.class_ = par.reverse_reaction_class(self.class_)
        self.forward_ts_graph, self.backward_ts_graph = (
            self.backward_ts_graph, self.forward_ts_graph)
        self.reactants_keys, self.products_keys = (
            self.products_keys, self.reactants_keys)

    def standardize_keys_(self):
        """ standardize keys and, optionally, sort reactant and product
        geometries in the standard order
        """
        rct_keys = list(map(sorted, self.reactants_keys))
        prd_keys = list(map(sorted, self.products_keys))
        rct_key_dct = {k: i for i, k in enumerate(itertools.chain(*rct_keys))}
        prd_key_dct = {k: i for i, k in enumerate(itertools.chain(*prd_keys))}
        self.reactants_keys = tuple(tuple(map(rct_key_dct.__getitem__, keys))
                                    for keys in self.reactants_keys)
        self.products_keys = tuple(tuple(map(prd_key_dct.__getitem__, keys))
                                   for keys in self.products_keys)
        self.forward_ts_graph = automol.graph.relabel(self.forward_ts_graph,
                                                      rct_key_dct)
        self.backward_ts_graph = automol.graph.relabel(self.backward_ts_graph,
                                                       prd_key_dct)

    def standardize_keys_and_sort_geometries_(self, rct_geos, prd_geos):
        """ standardize keys and line up geometries to match
        """
        rct_idxs, prd_idxs = self.sort_order()
        rct_geos = tuple(map(rct_geos.__getitem__, rct_idxs))
        prd_geos = tuple(map(prd_geos.__getitem__, prd_idxs))
        self.standardize_keys_()
        return rct_geos, prd_geos

    def insert_dummy_atoms_(self, dummy_key_dct, product=False):
        """ insert dummy atoms into the reactants or products

        :param dummy_key_dct: dummy atom keys, by key of the atom they are
            connected to
        :param product: insert dummy atoms into the products instead of the
            reactants?
        """
        if product:
            tsg = self.backward_ts_graph
            keys_lst = self.products_keys
        else:
            tsg = self.forward_ts_graph
            keys_lst = self.reactants_keys

        tsg = add_bonded_dummy_atoms(tsg, dummy_key_dct)
        keys_lst = list(map(list, keys_lst))
        for key, dummy_key in dummy_key_dct.items():
            for keys in keys_lst:
                if key in keys:
                    keys.append(dummy_key)
                    break
        keys_lst = tuple(map(tuple, map(sorted, keys_lst)))

        if product:
            self.backward_ts_graph = tsg
            self.products_keys = keys_lst
        else:
            self.forward_ts_graph = tsg
            self.reactants_keys = keys_lst

    def copy(self):
        """ return a copy of this Reaction
        """
        return Reaction(
            self.class_, self.forward_ts_graph, self.backward_ts_graph,
            self.reactants_keys, self.products_keys)

    def reactant_graphs(self):
        """ reactant graphs for this reaction, in order
        """
        rcts_gra = ts.reactants_graph(self.forward_ts_graph)
        rct_gras = [subgraph(rcts_gra, keys) for keys in self.reactants_keys]
        return tuple(rct_gras)

    def product_graphs(self):
        """ product graphs for this reaction, in order
        """
        prds_gra = ts.products_graph(self.forward_ts_graph)
        prd_gras = [subgraph(prds_gra, keys) for keys in self.products_keys]
        return tuple(prd_gras)

    def is_standardized(self):
        """ has this Reaction been standardized?
        """
        other = self.copy()
        other.standardize_keys_()
        return self == other

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
