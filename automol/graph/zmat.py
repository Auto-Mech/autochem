""" graph-based z-matrix builder
"""
from automol.graph._graph import longest_chain


def main(gra):
    """ v-matrix for a chain of heavy atoms
    """
    chain_keys = longest_chain(gra)
    print(chain_keys)


if __name__ == '__main__':
    import automol
    ICH = automol.smiles.inchi('CC')
    GEO = automol.inchi.geometry(ICH)
    GRA = automol.geom.connectivity_graph(GEO)
    main(GRA)
    # ZMA = automol.geom.zmatrix(GEO)
    # print(automol.zmatrix.string(ZMA, one_indexed=False))
    # print(automol.geom.zmatrix_torsion_coordinate_names(GEO))
