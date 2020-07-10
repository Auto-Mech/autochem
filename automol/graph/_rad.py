""" radical library
"""

from automol.graph._graph import unsaturated_atom_keys
from automol.graph._graph import add_atom_explicit_hydrogen_keys
from automol.graph._graph import remove_atoms
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import full_isomorphism
from automol.graph._graph_base import atom_symbols


def isomorphic_radical_graphs(gra):
    """ Generate a set of graphs that are isomorphic to a graph
        of a radical species
    """

    # Determine useful keys
    symbols = atom_symbols(gra)
    unsat_keys = unsaturated_atom_keys(gra)
    unsat_key = next(iter(unsat_keys))
    h_atm_key = max(symbols.keys()) + 1

    iso_gras = []
    for aidx, symbol in enumerate(symbols.values()):
        # Loop over saturated (non-radical) heavy atoms
        if symbol != 'H' and aidx != unsat_key:

            # Add hydrogen atom to radical atom
            new_graph = add_atom_explicit_hydrogen_keys(
                gra, {unsat_key: [h_atm_key]})

            # Remove hydrogen from saturated atom
            neighbors = atom_neighbor_keys(new_graph)
            for neigh in neighbors[aidx]:
                if symbols[neigh] == 'H':
                    aneighbor = neigh
                    break
            new_graph = remove_atoms(new_graph, [aneighbor])

            # Test to see if new radical species is the same as the original
            inv_atm_key_dct = full_isomorphism(gra, new_graph)
            if inv_atm_key_dct:
                iso_gras.append(new_graph)

    return iso_gras
