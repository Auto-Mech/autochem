""" construct transition state z-matrices
"""
import functools
import numpy
from qcelemental import constants as qcc
import automol.graph
import automol.graph.reaction
import automol.convert.zmatrix


def beta_scission(rct_zmas, prd_zmas):
    """ z-matrix for a beta-scission reaction
    """
    rct_zmas, rcts_gra = _standard_form_with_combined_graph(rct_zmas)
    prd_zmas, prds_gra = _standard_form_with_combined_graph(prd_zmas)
    rxn = automol.graph.reaction.beta_scission(rcts_gra, prds_gra)
    if rxn is not None:
        brk_bnd_key, = automol.graph.reaction.broken_bond_keys(rxn)
        ts_zma, = rct_zmas
        coo_dct = automol.zmatrix.coordinates(ts_zma)
        bnd_dist_coo_key = tuple(reversed(sorted(brk_bnd_key)))
        bnd_dist_name = next(coo_name for coo_name, coo_keys in coo_dct.items()
                             if bnd_dist_coo_key in coo_keys)
        ret = ts_zma, bnd_dist_name

    return ret


def addition(rct_zmas, prd_zmas):
    """ z-matrix for an addition reaction
    """
    bnd_dist_name = 'rts'
    bnd_dist_val = 3.

    rct_zmas, rcts_gra = _standard_form_with_combined_graph(rct_zmas)
    prd_zmas, prds_gra = _standard_form_with_combined_graph(prd_zmas)
    rxn = automol.graph.reaction.addition(rcts_gra, prds_gra)
    if rxn is not None:
        rct1_zma, rct2_zma = rct_zmas
        rct2_natms = automol.zmatrix.count(rct2_zma)

        frm_bnd_key, = automol.graph.reaction.formed_bond_keys(rxn)
        rct1_atm1_key, _ = sorted(frm_bnd_key)
        rct1_atm2_key, rct1_atm3_key, chain = _join_atom_keys(
            rct1_zma, rct1_atm1_key)

        if chain:
            join_val_dct = {
                bnd_dist_name: bnd_dist_val,
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
                bnd_dist_name: bnd_dist_val,
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
            [[bnd_dist_name, 'aabs1', 'babs1'],
             [None, 'aabs2', 'babs2'],
             [None, None, 'babs3']])[:rct2_natms]
        join_names[numpy.equal(join_keys, None)] = None

        join_name_set = set(numpy.ravel(join_names)) - {None}
        join_val_dct = {name: join_val_dct[name] for name in join_name_set}

        ts_zma = automol.zmatrix.join(
            rct1_zma, rct2_zma, join_keys, join_names, join_val_dct)
        ret = ts_zma, bnd_dist_name

    return ret


def hydrogen_abstraction(rct_zmas, prd_zmas):
    """ z-matrix for an addition reaction
    """
    bnd_dist_name = 'rts'
    bnd_dist_val = 3.
    rct_zmas, rcts_gra = _standard_form_with_combined_graph(rct_zmas)
    prd_zmas, prds_gra = _standard_form_with_combined_graph(prd_zmas)
    rxn = automol.graph.reaction.hydrogen_abstraction(rcts_gra, prds_gra)
    if rxn is not None:
        rct1_zma, rct2_zma = rct_zmas
        rct1_natms = automol.zmatrix.count(rct1_zma)
        rct2_natms = automol.zmatrix.count(rct2_zma)

        frm_bnd_key, = automol.graph.reaction.formed_bond_keys(rxn)
        brk_bnd_key, = automol.graph.reaction.broken_bond_keys(rxn)
        rct1_atm1_key = next(iter(frm_bnd_key & brk_bnd_key))
        rct1_atm2_key, rct1_atm3_key, chain = _join_atom_keys(
            rct1_zma, rct1_atm1_key)

        x_zma = ((('X', (None, None, None), (None, None, None)),), {})

        x_join_val_dct = {
            'rx': 1. * qcc.conversion_factor('angstrom', 'bohr'),
            'ax': 90. * qcc.conversion_factor('degree', 'radian'),
            'dx': 180. * qcc.conversion_factor('degree', 'radian'),
        }

        x_join_keys = numpy.array(
            [[rct1_atm1_key, rct1_atm2_key, rct1_atm3_key]])
        x_join_names = numpy.array([['rx', 'ax', 'dx']])
        x_join_names[numpy.equal(x_join_keys, None)] = None
        rct1_x_zma = automol.zmatrix.join(
            rct1_zma, x_zma, x_join_keys, x_join_names, x_join_val_dct)

        x_atm_key = rct1_natms

        if chain:
            join_val_dct = {
                bnd_dist_name: bnd_dist_val,
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
                bnd_dist_name: bnd_dist_val,
                'aabs1': 85. * qcc.conversion_factor('degree', 'radian'),
                'aabs2': 85. * qcc.conversion_factor('degree', 'radian'),
                'babs1': 180. * qcc.conversion_factor('degree', 'radian'),
                'babs2': 90. * qcc.conversion_factor('degree', 'radian'),
                'babs3': 90. * qcc.conversion_factor('degree', 'radian'),
            }

        join_keys = numpy.array(
            [[rct1_atm1_key, x_atm_key, rct1_atm2_key],
             [None, rct1_atm1_key, x_atm_key],
             [None, None, rct1_atm1_key]])[:rct2_natms]
        join_names = numpy.array(
            [[bnd_dist_name, 'aabs1', 'babs1'],
             [None, 'aabs2', 'babs2'],
             [None, None, 'babs3']])[:rct2_natms]
        join_names[numpy.equal(join_keys, None)] = None

        join_name_set = set(numpy.ravel(join_names)) - {None}
        join_val_dct = {name: join_val_dct[name] for name in join_name_set}

        ts_zma = automol.zmatrix.join(
            rct1_x_zma, rct2_zma, join_keys, join_names, join_val_dct)
        ret = ts_zma, bnd_dist_name

    return ret


def _standard_form_with_combined_graph(zmas):
    zmas = list(map(automol.zmatrix.standard_form, zmas))
    gras = list(map(automol.convert.zmatrix.graph, zmas))
    shift = 0
    for idx, (zma, gra) in enumerate(zip(zmas, gras)):
        zmas[idx] = automol.zmatrix.standard_form(zma, shift=shift)
        gras[idx] = automol.graph.transform_keys(gra, lambda x: x+shift)
        shift += len(automol.graph.atoms(gra))
    gra = functools.reduce(automol.graph.union, gras)
    gra = automol.graph.without_dummy_atoms(gra)
    zmas = tuple(zmas)
    return zmas, gra


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
