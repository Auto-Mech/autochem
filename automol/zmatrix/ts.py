""" construct transition state z-matrices
"""
import functools
import numpy
from qcelemental import constants as qcc
import automol.formula
import automol.graph
import automol.graph.trans
import automol.convert.zmatrix
import automol.zmatrix


def beta_scission(rct_zmas, prd_zmas):
    """ z-matrix for a beta-scission reaction
    """
    ret = None
    rct_zmas, rct_gras = _shifted_standard_forms_with_gaphs(rct_zmas)
    prd_zmas, prd_gras = _shifted_standard_forms_with_gaphs(prd_zmas)
    rcts_gra = functools.reduce(automol.graph.union, rct_gras)
    prds_gra = functools.reduce(automol.graph.union, prd_gras)
    tra = automol.graph.trans.beta_scission(rcts_gra, prds_gra)
    if tra is not None:
        brk_bnd_key, = automol.graph.trans.broken_bond_keys(tra)
        ts_zma, = rct_zmas
        coo_dct = automol.zmatrix.coordinates(ts_zma)
        dist_coo_key = tuple(reversed(sorted(brk_bnd_key)))
        dist_name = next(coo_name for coo_name, coo_keys in coo_dct.items()
                         if dist_coo_key in coo_keys)

        ts_name_dct = automol.zmatrix.standard_names(ts_zma)
        dist_name = ts_name_dct[dist_name]
        ts_zma = automol.zmatrix.standard_form(ts_zma)
        tors_names = automol.zmatrix.torsion_coordinate_names(ts_zma)

        ret = ts_zma, dist_name, tors_names

    return ret


def addition(rct_zmas, prd_zmas):
    """ z-matrix for an addition reaction
    """
    ret = None
    dist_name = 'rts'
    dist_val = 3.

    rct_zmas, rct_gras = _shifted_standard_forms_with_gaphs(rct_zmas)
    prd_zmas, prd_gras = _shifted_standard_forms_with_gaphs(prd_zmas)
    rcts_gra = functools.reduce(automol.graph.union, rct_gras)
    prds_gra = functools.reduce(automol.graph.union, prd_gras)
    tra = automol.graph.trans.addition(rcts_gra, prds_gra)
    if tra is not None:
        rct1_zma, rct2_zma = rct_zmas
        rct2_natms = automol.zmatrix.count(rct2_zma)

        frm_bnd_key, = automol.graph.trans.formed_bond_keys(tra)
        rct1_atm1_key, _ = sorted(frm_bnd_key)
        rct1_atm2_key, rct1_atm3_key, chain = _join_atom_keys(
            rct1_zma, rct1_atm1_key)

        join_val_dct = {
            dist_name: dist_val,
            'aabs1': 85. * qcc.conversion_factor('degree', 'radian'),
            'aabs2': 85. * qcc.conversion_factor('degree', 'radian'),
            'babs1': 90. * qcc.conversion_factor('degree', 'radian'),
            'babs2': 90. * qcc.conversion_factor('degree', 'radian'),
            'babs3': 90. * qcc.conversion_factor('degree', 'radian'),
        }

        join_keys = numpy.array(
            [[rct1_atm1_key, rct1_atm2_key, rct1_atm3_key],
             [None, rct1_atm1_key, rct1_atm2_key],
             [None, None, rct1_atm1_key]])[:rct2_natms]
        join_names = numpy.array(
            [[dist_name, 'aabs1', 'babs1'],
             [None, 'aabs2', 'babs2'],
             [None, None, 'babs3']])[:rct2_natms]
        join_names[numpy.equal(join_keys, None)] = None

        join_name_set = set(numpy.ravel(join_names)) - {None}
        join_val_dct = {name: join_val_dct[name] for name in join_name_set}

        ts_zma = automol.zmatrix.join(
            rct1_zma, rct2_zma, join_keys, join_names, join_val_dct)

        ts_name_dct = automol.zmatrix.standard_names(ts_zma)
        dist_name = ts_name_dct[dist_name]
        ts_zma = automol.zmatrix.standard_form(ts_zma)
        rct1_tors_names = automol.zmatrix.torsion_coordinate_names(rct1_zma)
        rct2_tors_names = automol.zmatrix.torsion_coordinate_names(rct2_zma)
        tors_names = (
            tuple(map(ts_name_dct.__getitem__, rct1_tors_names)) +
            tuple(map(ts_name_dct.__getitem__, rct2_tors_names))
        )

        if 'babs2' in ts_name_dct:
            tors_name = ts_name_dct['babs2']
            tors_names += (tors_name,)

        if 'babs3' in ts_name_dct:
            tors_name = ts_name_dct['babs3']
            tors_names += (tors_name,)

        ret = ts_zma, dist_name, tors_names

    return ret


def hydrogen_abstraction(rct_zmas, prd_zmas, sigma=False):
    """ z-matrix for an abstraction reaction
    """
    if sigma:
        ret = _sigma_hydrogen_abstraction(rct_zmas, prd_zmas)
    else:
        ret = _hydrogen_abstraction(rct_zmas, prd_zmas)
    return ret

def _sigma_hydrogen_abstraction(rct_zmas, prd_zmas):
    ret = None
    dist_name = 'rts'
    dist_val = 3.

    rxn_idxs = automol.formula.reac.argsort_hydrogen_abstraction(
        list(map(automol.convert.zmatrix.formula, rct_zmas)),
        list(map(automol.convert.zmatrix.formula, prd_zmas)))
    if rxn_idxs is not None:
        rct_idxs, prd_idxs = rxn_idxs
        rct_zmas = list(map(rct_zmas.__getitem__, rct_idxs))
        prd_zmas = list(map(prd_zmas.__getitem__, prd_idxs))
        rct_zmas, rct_gras = _shifted_standard_forms_with_gaphs(rct_zmas)
        prd_zmas, prd_gras = _shifted_standard_forms_with_gaphs(prd_zmas)
        rcts_gra = functools.reduce(automol.graph.union, rct_gras)
        prds_gra = functools.reduce(automol.graph.union, prd_gras)
        tra = automol.graph.trans.hydrogen_abstraction(rcts_gra, prds_gra)
        if tra is not None:
            rct1_gra, rct2_gra = rct_gras
            rct1_zma, rct2_zma = rct_zmas
            rct1_natms = automol.zmatrix.count(rct1_zma)
            rct2_natms = automol.zmatrix.count(rct2_zma)

            frm_bnd_key, = automol.graph.trans.formed_bond_keys(tra)
            brk_bnd_key, = automol.graph.trans.broken_bond_keys(tra)
            rct1_atm1_key = next(iter(frm_bnd_key & brk_bnd_key))
            rct2_atm1_key = next(iter(frm_bnd_key - brk_bnd_key))

            # if rct1 and rct2 are isomorphic, we may get an atom key on rct2.
            # in that case, determine the equivalent atom from rct1
            if rct1_atm1_key in automol.graph.atom_keys(rct2_gra):
                atm_key_dct = automol.graph.full_isomorphism(rct2_gra,
                                                             rct1_gra)
                assert atm_key_dct
                rct1_atm1_key = atm_key_dct[rct1_atm1_key]

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
            x_join_names = numpy.array([['rx', 'ax', 'dx']],
                                       dtype=numpy.object_)
            x_join_names[numpy.equal(x_join_keys, None)] = None
            x_join_name_set = set(numpy.ravel(x_join_names)) - {None}
            x_join_val_dct = {name: x_join_val_dct[name]
                              for name in x_join_name_set}
            rct1_x_zma = automol.zmatrix.join(
                rct1_zma, x_zma, x_join_keys, x_join_names, x_join_val_dct)

            rct1_x_atm_key = rct1_natms


            if rct2_atm1_key in automol.graph.atom_keys(rct1_gra):
                atm_key_dct = automol.graph.full_isomorphism(rct1_gra,
                                                             rct2_gra)
                assert atm_key_dct
                rct2_atm1_key = atm_key_dct[rct2_atm1_key]

            rct2_atm1_key -= rct1_natms
            print(automol.zmatrix.string(rct1_zma))
            print(automol.zmatrix.string(rct2_zma))
            print(rct2_atm1_key)
            assert rct2_atm1_key == 0
            # insert dummy atom as the second atom in reactant 2

            insert_keys = numpy.array(
                [[0, None, None],
                 [None, 1, None],
                 [None, None, 2]])[:rct2_natms]
            insert_names = numpy.array(
                [['rx2', None, None],
                 [None, 'ax2', None],
                 [None, None, 'dx2']])[:rct2_natms]
            insert_val_dct = {
                'rx2': 1. * qcc.conversion_factor('angstrom', 'bohr'),
                'ax2': 90. * qcc.conversion_factor('degree', 'radian'),
                'dx2': 180. * qcc.conversion_factor('degree', 'radian'),
            }
            insert_name_set = set(numpy.ravel(insert_names)) - {None}
            insert_val_dct = {name: insert_val_dct[name] for name in insert_name_set}
            rct2_x_zma = automol.zmatrix.insert_dummy_atom(
                rct2_zma, 1, insert_keys, insert_names, insert_val_dct
            )

            join_val_dct = {
                dist_name: dist_val,
                'aabs1': 85. * qcc.conversion_factor('degree', 'radian'),
                'aabs2': 85. * qcc.conversion_factor('degree', 'radian'),
                'babs1': 180. * qcc.conversion_factor('degree', 'radian'),
                'babs2': 90. * qcc.conversion_factor('degree', 'radian'),
                'babs3': 180. * qcc.conversion_factor('degree', 'radian'),
            }

            join_keys = numpy.array(
                [[rct1_atm1_key, rct1_x_atm_key, rct1_atm2_key],
                 [None, rct1_atm1_key, rct1_x_atm_key],
                 [None, None, rct1_atm1_key]])[:rct2_natms+1]
            join_names = numpy.array(
                [[dist_name, 'aabs1', 'babs1'],
                 [None, 'aabs2', 'babs2'],
                 [None, None, 'babs3']])[:rct2_natms+1]
            join_names[numpy.equal(join_keys, None)] = None

            join_name_set = set(numpy.ravel(join_names)) - {None}
            join_val_dct = {name: join_val_dct[name] for name in join_name_set}

            ts_zma = automol.zmatrix.join(
                rct1_x_zma, rct2_x_zma, join_keys, join_names, join_val_dct)

            ts_name_dct = automol.zmatrix.standard_names(ts_zma)
            print('babs test')
            print(ts_name_dct['babs2'])
            print(ts_name_dct['babs3'])
            inv_ts_name_dct = dict(map(reversed, ts_name_dct.items()))
            print(inv_ts_name_dct['D5'])
            ts_val_dct = automol.zmatrix.values(ts_zma)
            print(ts_val_dct['babs3'])
            dist_name = ts_name_dct[dist_name]
            ts_zma = automol.zmatrix.standard_form(ts_zma)
            rct1_tors_names = automol.zmatrix.torsion_coordinate_names(
                rct1_zma)
            rct2_tors_names = automol.zmatrix.torsion_coordinate_names(
                rct2_zma)
            tors_names = (
                tuple(map(ts_name_dct.__getitem__, rct1_tors_names)) +
                tuple(map(ts_name_dct.__getitem__, rct2_tors_names))
            )

            if 'babs2' in ts_name_dct:
                tors_name = ts_name_dct['babs2']
                tors_names += (tors_name,)

            if 'babs3' in ts_name_dct:
                tors_name = ts_name_dct['babs3']
                tors_names += (tors_name,)

            ret = ts_zma, dist_name, tors_names

    return ret


def _hydrogen_abstraction(rct_zmas, prd_zmas):
    ret = None
    dist_name = 'rts'
    dist_val = 3.

    rxn_idxs = automol.formula.reac.argsort_hydrogen_abstraction(
        list(map(automol.convert.zmatrix.formula, rct_zmas)),
        list(map(automol.convert.zmatrix.formula, prd_zmas)))
    if rxn_idxs is not None:
        rct_idxs, prd_idxs = rxn_idxs
        rct_zmas = list(map(rct_zmas.__getitem__, rct_idxs))
        prd_zmas = list(map(prd_zmas.__getitem__, prd_idxs))
        rct_zmas, rct_gras = _shifted_standard_forms_with_gaphs(rct_zmas)
        prd_zmas, prd_gras = _shifted_standard_forms_with_gaphs(prd_zmas)
        rcts_gra = functools.reduce(automol.graph.union, rct_gras)
        prds_gra = functools.reduce(automol.graph.union, prd_gras)
        tra = automol.graph.trans.hydrogen_abstraction(rcts_gra, prds_gra)
        if tra is not None:
            rct1_gra, rct2_gra = rct_gras
            rct1_zma, rct2_zma = rct_zmas
            rct1_natms = automol.zmatrix.count(rct1_zma)
            rct2_natms = automol.zmatrix.count(rct2_zma)

            frm_bnd_key, = automol.graph.trans.formed_bond_keys(tra)
            brk_bnd_key, = automol.graph.trans.broken_bond_keys(tra)
            rct1_atm1_key = next(iter(frm_bnd_key & brk_bnd_key))

            # if rct1 and rct2 are isomorphic, we may get an atom key on rct2.
            # in that case, determine the equivalent atom from rct1
            if rct1_atm1_key in automol.graph.atom_keys(rct2_gra):
                atm_key_dct = automol.graph.full_isomorphism(rct2_gra,
                                                             rct1_gra)
                assert atm_key_dct
                rct1_atm1_key = atm_key_dct[rct1_atm1_key]

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
            x_join_names = numpy.array([['rx', 'ax', 'dx']],
                                       dtype=numpy.object_)
            x_join_names[numpy.equal(x_join_keys, None)] = None
            x_join_name_set = set(numpy.ravel(x_join_names)) - {None}
            x_join_val_dct = {name: x_join_val_dct[name]
                              for name in x_join_name_set}
            rct1_x_zma = automol.zmatrix.join(
                rct1_zma, x_zma, x_join_keys, x_join_names, x_join_val_dct)

            x_atm_key = rct1_natms

            join_val_dct = {
                dist_name: dist_val,
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
                [[dist_name, 'aabs1', 'babs1'],
                 [None, 'aabs2', 'babs2'],
                 [None, None, 'babs3']])[:rct2_natms]
            join_names[numpy.equal(join_keys, None)] = None

            join_name_set = set(numpy.ravel(join_names)) - {None}
            join_val_dct = {name: join_val_dct[name] for name in join_name_set}

            ts_zma = automol.zmatrix.join(
                rct1_x_zma, rct2_zma, join_keys, join_names, join_val_dct)

            ts_name_dct = automol.zmatrix.standard_names(ts_zma)
            dist_name = ts_name_dct[dist_name]
            ts_zma = automol.zmatrix.standard_form(ts_zma)
            rct1_tors_names = automol.zmatrix.torsion_coordinate_names(
                rct1_zma)
            rct2_tors_names = automol.zmatrix.torsion_coordinate_names(
                rct2_zma)
            tors_names = (
                tuple(map(ts_name_dct.__getitem__, rct1_tors_names)) +
                tuple(map(ts_name_dct.__getitem__, rct2_tors_names))
            )

            if 'babs2' in ts_name_dct:
                tors_name = ts_name_dct['babs2']
                tors_names += (tors_name,)

            if 'babs3' in ts_name_dct:
                tors_name = ts_name_dct['babs3']
                tors_names += (tors_name,)

            ret = ts_zma, dist_name, tors_names

    return ret


def _shifted_standard_forms_with_gaphs(zmas):
    gras = list(map(automol.convert.zmatrix.graph, zmas))
    shift = 0
    for idx, (zma, gra) in enumerate(zip(zmas, gras)):
        zmas[idx] = automol.zmatrix.standard_form(zma, shift=shift)
        gras[idx] = automol.graph.transform_keys(gra, lambda x: x+shift)
        shift += len(automol.graph.atoms(gra))
    zmas = tuple(zmas)
    gras = tuple(map(automol.graph.without_dummy_atoms, gras))
    return zmas, gras


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
