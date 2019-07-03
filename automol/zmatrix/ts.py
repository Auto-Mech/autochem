""" construct transition state z-matrices
"""
import functools
import automol.graph
import automol.graph.reaction
import automol.convert.zmatrix


# def beta_scission(rct_zmas, prd_zmas):
#     """ z-matrix for a beta scission reaction
#     """


def addition(rct_zmas, prd_zmas):
    """ z-matrix for an addition reaction
    """
    rcts_gra = _combined_graph(rct_zmas)
    prds_gra = _combined_graph(prd_zmas)
    rxn = automol.graph.reaction.addition(rcts_gra, prds_gra)
    frm_bnd_key, = automol.graph.reaction.formed_bond_keys(rxn)
    print(rcts_gra)
    print(prds_gra)
    print(frm_bnd_key)


def _combined_graph(zmas):
    gras = list(map(automol.convert.zmatrix.graph, zmas))
    gras = list(map(automol.graph.without_dummy_atoms, gras))
    shift = 0
    for idx, gra in enumerate(gras):
        gras[idx] = automol.graph.transform_keys(gra, lambda x: x+shift)
        shift += len(automol.graph.atoms(gra))
    gra = functools.reduce(automol.graph.union, gras)
    return gra


if __name__ == '__main__':
    RCT_ZMAS = [
        ((('H', (None, None, None), (None, None, None)),),
         {}),
        ((('O', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None))),
         {'R1': 2.15608}),
    ]
    PRD_ZMAS = [
        ((('O', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None))),
         {'R1': 2.48959, 'R2': 1.86213, 'A2': 1.9084302705931997})
    ]
    addition(RCT_ZMAS, PRD_ZMAS)
    # addition(reversed(RCT_ZMAS), PRD_ZMAS)
