""" construct transition state z-matrices
"""
import functools
import numpy
from qcelemental import constants as qcc
import automol.graph
import automol.graph.reaction
import automol.convert.zmatrix


# def beta_scission(rct_zmas, prd_zmas):
#     """ z-matrix for a beta scission reaction
#     """


def _join_atom_keys(zma, atm1_key):
    """ returns available join atom keys (if available) and a boolean
    indicating whether the atoms are in a chain or not
    """
    gra = automol.convert.zmatrix.graph(zma)
    atm1_chain = (
        automol.graph.atom_longest_chains(gra)[atm1_key])
    atm1_ngb_keys = (
        automol.graph.atom_neighbor_keys(gra)[atm1_key])
    if len(atm1_chain) == 1:
        atm2_key = None
        atm3_key = None
        chain = False
    elif len(atm1_chain) == 2 and len(atm1_ngb_keys) == 1:
        atm2_key = atm1_chain[1]
        atm3_key = None
        chain = False
    elif len(atm1_chain) == 2:
        atm2_key = atm1_chain[1]
        atm3_key = sorted(atm1_ngb_keys - {atm2_key})[0]
        chain = False
    else:
        atm2_key, atm3_key = atm1_chain[1:3]
        chain = True

    return atm2_key, atm3_key, chain


def addition(rct_zmas, prd_zmas):
    """ z-matrix for an addition reaction
    """
    rcts_gra = _combined_graph(rct_zmas)
    prds_gra = _combined_graph(prd_zmas)
    rxn = automol.graph.reaction.addition(rcts_gra, prds_gra)
    if rxn is not None:
        frm_bnd_key, = automol.graph.reaction.formed_bond_keys(rxn)
        rct1_atm1_key, _ = sorted(frm_bnd_key)
        rct1_zma, rct2_zma = rct_zmas

        rct2_natms = automol.zmatrix.count(rct2_zma)

        rct1_atm2_key, rct1_atm3_key, chain = _join_atom_keys(
            rct1_zma, rct1_atm1_key)
        print(automol.zmatrix.string(rct1_zma))
        print(automol.zmatrix.string(rct2_zma))

        if chain:
            join_val_dct = {
                'rts': 5.,
                'aabs1': 85. * qcc.conversion_factor('degree', 'radian'),
                'aabs2': 85. * qcc.conversion_factor('degree', 'radian'),
                'babs1': 180. * qcc.conversion_factor('degree', 'radian'),
                'babs2': 90. * qcc.conversion_factor('degree', 'radian'),
                'babs3': 90. * qcc.conversion_factor('degree', 'radian'),
            }
        else:
            # the same for now -- eventually, we should have different defaults
            # in this case
            join_val_dct = {
                'rts': 5.,
                'aabs1': 85. * qcc.conversion_factor('degree', 'radian'),
                'aabs2': 85. * qcc.conversion_factor('degree', 'radian'),
                'babs1': 180. * qcc.conversion_factor('degree', 'radian'),
                'babs2': 90. * qcc.conversion_factor('degree', 'radian'),
                'babs3': 90. * qcc.conversion_factor('degree', 'radian'),
            }

        join_keys = numpy.array(
            [[rct1_atm1_key, rct1_atm2_key, rct1_atm3_key],
             [None, rct1_atm1_key, rct1_atm2_key],
             [None, None, rct1_atm1_key]])[:rct2_natms]
        join_names = numpy.array(
            [['rts', 'aabs1', 'babs1'],
             [None, 'aabs2', 'babs2'],
             [None, None, 'babs3']])[:rct2_natms]
        join_names[numpy.equal(join_keys, None)] = None

        join_name_set = set(numpy.ravel(join_names)) - {None}
        join_val_dct = {name: join_val_dct[name] for name in join_name_set}

        print(numpy.array(join_keys))
        print(numpy.array(join_names))

        ts_zma = automol.zmatrix.join(
            rct1_zma, rct2_zma, join_keys, join_names, join_val_dct)

        print(automol.zmatrix.string(ts_zma))


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
        ((('O', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None))),
         {'R1': 2.15608}),
        ((('H', (None, None, None), (None, None, None)),),
         {}),
    ]
    PRD_ZMAS = [
        ((('O', (None, None, None), (None, None, None)),
          ('O', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None))),
         {'R1': 2.48959, 'R2': 1.86213, 'A2': 1.9084302705931997})
    ]
    addition(RCT_ZMAS, PRD_ZMAS)
    # addition(reversed(RCT_ZMAS), PRD_ZMAS)
