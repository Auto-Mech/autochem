""" graph geometry library

For generating heuristic coordinates and z-matrices from graphs
"""
import automol.zmat
from automol.graph._graph import atom_symbols
from automol.graph._graph import implicit
from automol.graph._graph import longest_chain


def main_function(gra):
    """ main function that does stuff
    """
    gra = implicit(gra)

    sym_dct = atom_symbols(gra)
    chain = longest_chain(gra)

    # start from an empty z-matrix and build up
    zma = automol.zmat.add_atom((), sym_dct[chain[0]], (), ())

    print(automol.zmat.string(zma))
    print(chain)


if __name__ == '__main__':
    import automol

    ICH = automol.smiles.inchi('CCCC')
    GRA = automol.inchi.graph(ICH)
    print(automol.graph.string(GRA, one_indexed=False))
    main_function(GRA)
