""" utilities
"""
import automol


def standardize_reagent_graphs_keys(gras):
    """ standardize keys for a series of reagent graphs
    """
    ret_gras = []
    shift = 0
    for gra in gras:
        gra = automol.graph.standard_keys(gra)
        gra = automol.graph.transform_keys(gra, lambda x: x + shift)
        ret_gras.append(gra)
        shift += automol.formula.atom_count(automol.graph.formula(gra))
    return tuple(ret_gras)
