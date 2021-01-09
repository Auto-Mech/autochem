""" graph-based reaction classification

Function arguments:
    Each function takes a list of reactant graphs and a list of product graphs.
    Note that the reactant graphs *cannot* have overlapping atom keys, and
    likewise for the product graphs. Otherwise, there would be no way to
    express the bonds broken and formed between reactants.
"""
import itertools
import numpy
import automol.convert.graph
from automol import par
from automol.graph import ts
from automol.graph._graph_base import atom_keys
from automol.graph._graph_base import string
from automol.graph._graph import explicit
from automol.graph._graph import full_isomorphism
from automol.graph._graph import union_from_sequence
from automol.graph._graph import unsaturated_atom_keys
from automol.graph._graph import without_stereo_parities
from automol.graph._graph import add_atom_explicit_hydrogen_keys


class GraphReaction:
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
        """ contructor
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
        assert full_isomorphism(ts.reverse(forw_tsg), back_tsg)

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
        rct_keys = list(map(tuple, map(sorted, self.reactants_keys)))
        prd_keys = list(map(tuple, map(sorted, self.products_keys)))
        rct_idxs = tuple(numpy.argsort(numpy.array(rct_keys, dtype=object)))
        prd_idxs = tuple(numpy.argsort(numpy.array(prd_keys, dtype=object)))
        return rct_idxs, prd_idxs

    def standardize_keys(self):
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


def trivial(rct_gras, prd_gras):
    """ find a trivial reaction, with the same reactants and products
    """
    _assert_is_valid_reagent_graph_list(rct_gras)
    _assert_is_valid_reagent_graph_list(prd_gras)

    rxns = []

    if len(rct_gras) == len(prd_gras):
        prd_gras = list(prd_gras)

        rct_idxs = []
        prd_idxs = []

        # One at a time, find matches for each reactant; track the positions to
        # get the right sort order
        for rct_idx, rct_gra in enumerate(rct_gras):
            prd_idx = next((idx for idx, prd_gra in enumerate(prd_gras)
                            if full_isomorphism(rct_gra, prd_gra)), None)

            if prd_idx is not None:
                rct_idxs.append(rct_idx)
                prd_idxs.append(prd_idx)
                prd_gras.pop(prd_idx)
            else:
                break

        if rct_idxs and prd_idxs:
            # reorder the reactants and products
            rct_gras = list(map(rct_gras.__getitem__, rct_idxs))
            prd_gras = list(map(prd_gras.__getitem__, prd_idxs))

            rcts_gra = union_from_sequence(rct_gras)
            prds_gra = union_from_sequence(prd_gras)

            rxns.append(GraphReaction(
                rxn_cls=par.ReactionClass.TRIVIAL,
                forw_tsg=ts.graph(rcts_gra, [], []),
                back_tsg=ts.graph(prds_gra, [], []),
                rcts_keys=list(map(atom_keys, rct_gras)),
                prds_keys=list(map(atom_keys, prd_gras)),
            ))

    return tuple(rxns)


def hydrogen_abstractions(rct_gras, prd_gras):
    """ find hydrogen abstractions consistent with these reactants and products

    Hydrogen abstractions are identified first by checking whether the
    molecular formulas are consistent with a reaction of the form R1H + R2 =>
    R2H + R1. If they do, we identify the abstraction sites by adding hydrogens
    to unsaturated sites of the R1 product to see if we get the R1H reactant.
    We then do the same for the R2 reactant and the R2H product.
    """
    _assert_is_valid_reagent_graph_list(rct_gras)
    _assert_is_valid_reagent_graph_list(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 2:
        rct_fmls = list(map(automol.convert.graph.formula, rct_gras))
        prd_fmls = list(map(automol.convert.graph.formula, prd_gras))

        ret = automol.formula.reac.argsort_hydrogen_abstraction(
            rct_fmls, prd_fmls)
        if ret:
            rct_idxs_, prd_idxs_ = ret
            rct_gras = list(map(rct_gras.__getitem__, rct_idxs_))
            prd_gras = list(map(prd_gras.__getitem__, prd_idxs_))

            q1h_gra, q2_gra = rct_gras
            q2h_gra, q1_gra = prd_gras

            rets1 = _partial_hydrogen_abstraction(q1h_gra, q1_gra)
            rets2 = _partial_hydrogen_abstraction(q2h_gra, q2_gra)
            for ret1, ret2 in itertools.product(rets1, rets2):
                f_q1h_q_atm_key, f_q1h_h_atm_key, b_q2_q_atm_key = ret1
                b_q1h_q_atm_key, b_q1h_h_atm_key, f_q2_q_atm_key = ret2

                # Create the forward ts graph
                rcts_gra = union_from_sequence(rct_gras)
                f_frm_bnd_key = frozenset({f_q2_q_atm_key, f_q1h_h_atm_key})
                f_brk_bnd_key = frozenset({f_q1h_q_atm_key, f_q1h_h_atm_key})
                forw_tsg = ts.graph(rcts_gra,
                                    frm_bnd_keys=[f_frm_bnd_key],
                                    brk_bnd_keys=[f_brk_bnd_key])

                # Create the backward ts graph
                prds_gra = union_from_sequence(prd_gras)
                b_frm_bnd_key = frozenset({b_q2_q_atm_key, b_q1h_h_atm_key})
                b_brk_bnd_key = frozenset({b_q1h_q_atm_key, b_q1h_h_atm_key})
                back_tsg = ts.graph(prds_gra,
                                    frm_bnd_keys=[b_frm_bnd_key],
                                    brk_bnd_keys=[b_brk_bnd_key])

                # Create the reaction object
                rxns.append(GraphReaction(
                    rxn_cls=par.ReactionClass.HYDROGEN_ABSTRACTION,
                    forw_tsg=forw_tsg,
                    back_tsg=back_tsg,
                    rcts_keys=list(map(atom_keys, rct_gras)),
                    prds_keys=list(map(atom_keys, prd_gras)),
                ))

    return tuple(rxns)


def _partial_hydrogen_abstraction(qh_gra, q_gra):
    rets = []

    h_atm_key = max(atom_keys(q_gra)) + 1
    uns_atm_keys = unsaturated_atom_keys(q_gra)
    for atm_key in uns_atm_keys:
        q_gra_h = add_atom_explicit_hydrogen_keys(
            q_gra, {atm_key: [h_atm_key]})
        inv_atm_key_dct = full_isomorphism(q_gra_h, qh_gra)
        if inv_atm_key_dct:
            qh_q_atm_key = inv_atm_key_dct[atm_key]
            qh_h_atm_key = inv_atm_key_dct[h_atm_key]
            q_q_atm_key = atm_key
            rets.append((qh_q_atm_key, qh_h_atm_key, q_q_atm_key))

    return rets


def _assert_is_valid_reagent_graph_list(gras):
    gras_str = '\n---\n'.join(map(string, gras))
    assert _are_all_explicit(gras), (
        "Implicit hydrogens are not allowed here!\nGraphs:\n{}"
        .format(gras_str))
    assert _have_no_stereo_assignments(gras), (
        "Stereo assignments are not allowed here!\nGraphs:\n{}"
        .format(gras_str))
    assert _have_no_common_atom_keys(gras), (
        "Overlapping atom keys are not allowed here!\nGraphs:\n{}"
        .format(gras_str))


def _are_all_explicit(gras):
    return all(gra == explicit(gra) for gra in gras)


def _have_no_stereo_assignments(gras):
    return all(gra == without_stereo_parities(gra) for gra in gras)


def _have_no_common_atom_keys(gras):
    atm_keys = list(itertools.chain(*map(atom_keys, gras)))
    return len(atm_keys) == len(set(atm_keys))
