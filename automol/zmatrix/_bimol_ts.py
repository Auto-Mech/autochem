""" construct transition state z-matrices
"""

import functools
import numpy
from qcelemental import constants as qcc
import automol
from automol.zmatrix._util import shifted_standard_zmas_graphs
from automol.zmatrix._util import join_atom_keys
from automol.zmatrix._util import include_babs3
from automol.zmatrix._util import reorder_zma_for_radicals


ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')


def insertion(rct_zmas, prd_zmas):
    """ z-matrix for an insertion reaction
    """
    ret = None
    dist_name = 'rts'
    dist_val = 3.
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    tras, rct_idxs, _ = automol.graph.reac.insertion(rct_gras, prd_gras)
    if tras:
        tra = tras[0]
        frm_bnd_key, _ = automol.graph.trans.formed_bond_keys(tra)

        # get the atom on react 1 that is being attacked (bond is forming)
        rct1_atm1_key = list(frm_bnd_key)[0]

        # figure out atoms in the chain to define the dummy atom
        _, rct2_gra = map(rct_gras.__getitem__, rct_idxs)
        rct1_zma, rct2_zma = map(rct_zmas.__getitem__, rct_idxs)
        rct1_natms = automol.zmatrix.count(rct1_zma)
        rct2_natms = automol.zmatrix.count(rct2_zma)

        rct1_atm2_key, rct1_atm3_key, _ = join_atom_keys(
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
            'babs1': 170. * qcc.conversion_factor('degree', 'radian'),
            'babs2': 85. * qcc.conversion_factor('degree', 'radian'),
            'babs3': 85. * qcc.conversion_factor('degree', 'radian'),
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
        # CHECK IF THE INSERTION MESSES WITH THE TORSION ASSIGNMENTS
        rct1_tors_names = automol.zmatrix.torsion_coordinate_names(
            rct1_zma)
        rct2_tors_names = automol.zmatrix.torsion_coordinate_names(
            rct2_zma)
        tors_names = (
            tuple(map(ts_name_dct.__getitem__, rct1_tors_names)) +
            tuple(map(ts_name_dct.__getitem__, rct2_tors_names))
        )

        if 'babs2' in ts_name_dct:
            geo1 = automol.convert.zmatrix.geometry(rct1_zma)
            if not automol.geom.is_linear(geo1):
                tors_name = ts_name_dct['babs2']
                tors_names += (tors_name,)

        if 'babs3' in ts_name_dct and include_babs3(frm_bnd_key, rct2_gra):
            tors_name = ts_name_dct['babs3']
            tors_names += (tors_name,)
        ret = ts_zma, dist_name, tors_names

    return ret


def substitution(rct_zmas, prd_zmas):
    """ z-matrix for a substitution reaction
    Presume that radical substitutions at a pi bond occur instead as
    a sequence of addition and elimination.
    Also, for now presume that we are only interested in
    radical molecule substitutions
    """
    ret = None

    # Set the name and value for the bond being formed
    dist_name = 'rts'
    dist_val = 3.

    # first determine if this is a radical molecule reaction and then
    # if the radical is the first species
    # reorder to put it second
    rad_cnt = 0
    mol_cnt = 0
    for idx, rct_zma in enumerate(rct_zmas):
        rad_keys = automol.graph.resonance_dominant_radical_atom_keys(
            automol.geom.graph(automol.zmatrix.geometry(rct_zma)))
        ich = automol.geom.inchi(automol.zmatrix.geometry(rct_zma))
        is_co = (ich == 'InChI=1S/CO/c1-2')
        if rad_keys and not is_co:
            rad_idx = idx
            rad_cnt += 1
        else:
            # mol_idx = idx
            mol_cnt += 1

    if rad_cnt == 1 and mol_cnt == 1:
        if rad_idx == 0:
            rct2_zma, rct1_zma = rct_zmas
            rct_zmas = [rct1_zma, rct2_zma]
        # Confirm the reaction type and build the appropriate Z-Matrix

        rct2_gra = automol.zmatrix.graph(rct_zmas[1], remove_stereo=True)
        rad_atm_keys = automol.graph.resonance_dominant_radical_atom_keys(
            rct2_gra)
        if 0 not in rad_atm_keys:
            rct_zmas[1] = reorder_zma_for_radicals(
                rct_zmas[1], min(rad_atm_keys))
            rct2_gra = automol.zmatrix.graph(rct_zmas[1], remove_stereo=True)
            rad_atm_keys = automol.graph.resonance_dominant_radical_atom_keys(
                rct2_gra)
            # following assert checks to ensure that
            # the first atom in the second reactant is a radical
            # this is required for the remainder of the routine
            assert 0 in rad_atm_keys

        rct_zmas, rct_gras = shifted_standard_zmas_graphs(
            rct_zmas, remove_stereo=True)
        prd_zmas, prd_gras = shifted_standard_zmas_graphs(
            prd_zmas, remove_stereo=True)
        rcts_gra = functools.reduce(automol.graph.union, rct_gras)
        prds_gra = functools.reduce(automol.graph.union, prd_gras)

        tras, rct_idxs, _ = automol.graph.reac.substitution(rcts_gra, prds_gra)
    else:
        tras = []

    if tras:
        tra = tras[0]
        frm_bnd_key, = automol.graph.trans.formed_bond_keys(tra)

        # get the atom on react 1 that is being attacked (bond is forming)
        rct1_atm1_key = list(frm_bnd_key)[0]

        rct1_zma, rct2_zma = map(rct_zmas.__getitem__, rct_idxs)
        _, rct2_gra = map(rct_gras.__getitem__, rct_idxs)
        rct1_natms = automol.zmatrix.count(rct1_zma)
        rct2_natms = automol.zmatrix.count(rct2_zma)

        # figure out atoms in the chain to define the dummy atom
        rct1_atm2_key, rct1_atm3_key, _ = join_atom_keys(
            rct1_zma, rct1_atm1_key)

        # Join the reactant 1 ZMAT with the dummy atomx
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

        # Join the React 2 ZMAT with the reac1_x ZMAT
        x_atm_key = rct1_natms

        join_val_dct = {
            dist_name: dist_val,
            'aabs1': 85. * qcc.conversion_factor('degree', 'radian'),
            'aabs2': 85. * qcc.conversion_factor('degree', 'radian'),
            'babs1': 170. * qcc.conversion_factor('degree', 'radian'),
            'babs2': 85. * qcc.conversion_factor('degree', 'radian'),
            'babs3': 85. * qcc.conversion_factor('degree', 'radian'),
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

        # Get the names of the coordinates of the breaking and forming bond
        ts_name_dct = automol.zmatrix.standard_names(ts_zma)
        form_dist_name = ts_name_dct[dist_name]

        # Get the torsional coordinates of the transition state
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
            geo1 = automol.convert.zmatrix.geometry(rct1_zma)
            if not automol.geom.is_linear(geo1):
                tors_name = ts_name_dct['babs2']
                tors_names += (tors_name,)

        if 'babs3' in ts_name_dct and include_babs3(frm_bnd_key, rct2_gra):
            tors_name = ts_name_dct['babs3']
            tors_names += (tors_name,)

        # Set the reactant graph

        # Set info to be returned
        ret = ts_zma, form_dist_name, tors_names

    return ret


def addition(rct_zmas, prd_zmas, rct_tors=()):
    """ z-matrix for an addition reaction
    """
    ret = None
    dist_name = 'rts'
    dist_val = 3.

    count1 = automol.zmatrix.count(rct_zmas[0])
    if len(rct_zmas) == 2:
        count2 = automol.zmatrix.count(rct_zmas[1])
        if count1 == 1 or count1 < count2:
            rct2_zma, rct1_zma = rct_zmas
            rct_zmas = [rct1_zma, rct2_zma]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    tras, _, _ = automol.graph.reac.addition(rct_gras, prd_gras)
    # print('tras')
    # for tra in tras:
    if tras:
        tra = tras[0]
        rct1_zma, rct2_zma = rct_zmas
        rct1_gra, rct2_gra = rct_gras
        rct2_natms = automol.zmatrix.count(rct2_zma)

        frm_bnd_key, = automol.graph.trans.formed_bond_keys(tra)
        rct1_atm1_key, _ = sorted(frm_bnd_key)

        # Replace old routine with new one based on unsaturated atoms
        # Start of new routine
        rct1_unsat_keys = automol.graph.unsaturated_atom_keys(rct1_gra)
        # get neighbor keys for rct1_atm1_key
        neighbor_dct = automol.graph.atom_neighbor_keys(rct1_gra)
        atm1_nghbr_keys = neighbor_dct[rct1_atm1_key]
        # shift keys for unsaturated atoms and for atm1_nghbr keys
        # to include dummy atoms

        # second atom key
        # choose rct1_atm2 as first unsaturated neighbor to rct1_atm1
        rct1_atm2_key = None
        for key in atm1_nghbr_keys:
            if key in rct1_unsat_keys:
                rct1_atm2_key = key
                break

        # if no unsaturated neighbor, choose any atm1_nghbr
        if rct1_atm2_key is None:
            for key in atm1_nghbr_keys:
                rct1_atm2_key = key
                break

        # third atom key
        # if atm1 has two neighbors choose third atom key as second neighbor
        rct1_atm3_key = []
        for key in atm1_nghbr_keys:
            if key != rct1_atm2_key:
                rct1_atm3_key = key
                break

        # else choose it as second neighbor to atm 2
        if not rct1_atm3_key:
            atm2_nghbr_keys = neighbor_dct[rct1_atm2_key]
            for key in atm2_nghbr_keys:
                if key != rct1_atm1_key:
                    rct1_atm3_key = key
                    break

        # check if rct1_atm1, rct1_atm2, rct1_atm3 are colinear
        # if they are choose a dummy atom
        rct1_geo = automol.zmatrix.geometry(rct1_zma)
        if rct1_atm3_key:
            rct1_sub_geo = (
                rct1_geo[rct1_atm1_key], rct1_geo[rct1_atm2_key],
                rct1_geo[rct1_atm3_key])
        else:
            rct1_sub_geo = (rct1_geo[rct1_atm1_key], rct1_geo[rct1_atm2_key])
        if automol.geom.is_linear(rct1_sub_geo) or not rct1_atm3_key:
            # have to regenerate atm2_nghbr_keys to include dummy atom keys!
            rct1_key_mat = automol.zmatrix.key_matrix(rct1_zma)
            rct1_keys = [row[0] for row in rct1_key_mat]
            if rct1_keys[rct1_atm2_key] is not None:
                atm2_nghbr_keys = [rct1_keys[rct1_atm2_key]]
            else:
                atm2_nghbr_keys = []
            for idx, rct1_key in enumerate(rct1_keys):
                if rct1_key == rct1_atm2_key:
                    atm2_nghbr_keys.append(idx)
            new_atm3 = False
            for atm_key in atm2_nghbr_keys:
                if atm_key not in (rct1_atm3_key, rct1_atm1_key):
                    rct1_atm3_key = atm_key
                    new_atm3 = True
                    break
            if not new_atm3:
                # if there are no dummy atoms connected to the rct1_atm2 then
                # search for dummy atom keys connected to rct1_atm1
                # now regenerate atm1_nghbr_keys to include dummy atom keys!
                if rct1_keys[rct1_atm1_key] is not None:
                    atm1_nghbr_keys = [rct1_keys[rct1_atm1_key]]
                else:
                    atm1_nghbr_keys = []
                for idx, rct1_key in enumerate(rct1_keys):
                    if rct1_key == rct1_atm1_key:
                        atm1_nghbr_keys.append(idx)
                new_atm3 = False
                for atm_key in atm1_nghbr_keys:
                    if atm_key not in (rct1_atm3_key, rct1_atm2_key):
                        rct1_atm3_key = atm_key
                        new_atm3 = True
                        break

        join_val_dct = {
            dist_name: dist_val,
            'aabs1': 85. * qcc.conversion_factor('degree', 'radian'),
            'aabs2': 85. * qcc.conversion_factor('degree', 'radian'),
            'babs1': 85. * qcc.conversion_factor('degree', 'radian'),
            'babs2': 85. * qcc.conversion_factor('degree', 'radian'),
            'babs3': 85. * qcc.conversion_factor('degree', 'radian'),
        }

        rct1_natom = automol.zmatrix.count(rct1_zma)
        rct2_natom = automol.zmatrix.count(rct2_zma)

        if rct1_natom == 1 and rct2_natom == 1:
            raise NotImplementedError
        if rct1_natom == 2 and rct2_natom == 1:
            join_keys = numpy.array(
                [[rct1_atm1_key, rct1_atm2_key, None]])
            join_names = numpy.array(
                [[dist_name, 'aabs1', None]])
        elif rct1_natom == 2 and rct2_natom == 2:
            join_keys = numpy.array(
                [[rct1_atm1_key, rct1_atm2_key, rct1_atm3_key],
                 [None, rct1_atm1_key, rct1_atm2_key]])
            join_names = numpy.array(
                [[dist_name, 'aabs1', None],
                 [None, 'aabs2', 'babs2']])
        else:
            join_keys = numpy.array(
                [[rct1_atm1_key, rct1_atm2_key, rct1_atm3_key],
                 [None, rct1_atm1_key, rct1_atm2_key],
                 [None, None, rct1_atm1_key]])[:rct2_natms]
            join_names = numpy.array(
                [[dist_name, 'aabs1', 'babs1'],
                 [None, 'aabs2', 'babs2'],
                 [None, None, 'babs3']])[:rct2_natms]

        join_name_set = set(numpy.ravel(join_names)) - {None}
        join_val_dct = {name: join_val_dct[name] for name in join_name_set}

        ts_zma = automol.zmatrix.join(
            rct1_zma, rct2_zma, join_keys, join_names, join_val_dct)

        ts_name_dct = automol.zmatrix.standard_names(ts_zma)
        dist_name = ts_name_dct[dist_name]
        ts_zma = automol.zmatrix.standard_form(ts_zma)
        rct1_tors_names = automol.zmatrix.torsion_coordinate_names(rct1_zma)

        if rct_tors:
            rct2_tors_names = rct_tors
        else:
            rct2_tors_names = automol.zmatrix.torsion_coordinate_names(
                rct2_zma)
        tors_names = (
            tuple(map(ts_name_dct.__getitem__, rct1_tors_names)) +
            tuple(map(ts_name_dct.__getitem__, rct2_tors_names))
        )

        if 'babs2' in ts_name_dct:
            tors_name = ts_name_dct['babs2']
            tors_names += (tors_name,)

        if 'babs3' in ts_name_dct and include_babs3(frm_bnd_key, rct2_gra):
            tors_name = ts_name_dct['babs3']
            tors_names += (tors_name,)

        # print('frm_bnd_key test before shift:', frm_bnd_key)
        frm_bnd_key = shift_vals_from_dummy(frm_bnd_key, ts_zma)
        # print('frm_bnd_key test after shift:', frm_bnd_key)

        # Build reactants graph
        rcts_gra = automol.graph.union_from_sequence(rct_gras)

        ret = ts_zma, dist_name, frm_bnd_key, tors_names, rcts_gra

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

    # in the future, we can get the rxn_idxs directly from
    # graph.reac.hydrogen_abstraction, because it returns them now
    rxn_idxs = automol.formula.reac.argsort_hydrogen_abstraction(
        list(map(automol.convert.zmatrix.formula, rct_zmas)),
        list(map(automol.convert.zmatrix.formula, prd_zmas)))
    if rxn_idxs is not None:
        rct_idxs, prd_idxs = rxn_idxs
        rct_zmas = list(map(rct_zmas.__getitem__, rct_idxs))
        prd_zmas = list(map(prd_zmas.__getitem__, prd_idxs))
        rct_zmas, rct_gras = shifted_standard_zmas_graphs(
            rct_zmas, remove_stereo=True)
        prd_zmas, prd_gras = shifted_standard_zmas_graphs(
            prd_zmas, remove_stereo=True)
        tras, _, _ = automol.graph.reac.hydrogen_abstraction(
            rct_gras, prd_gras)
        if tras:
            tra = tras[0]
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

            rct1_atm2_key, rct1_atm3_key, _ = join_atom_keys(
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
            insert_val_dct = {name: insert_val_dct[name]
                              for name in insert_name_set}
            rct2_x_zma = automol.zmatrix.insert_dummy_atom(
                rct2_zma, 1, insert_keys, insert_names, insert_val_dct
            )

            join_val_dct = {
                dist_name: dist_val,
                'aabs1': 85. * qcc.conversion_factor('degree', 'radian'),
                'aabs2': 85. * qcc.conversion_factor('degree', 'radian'),
                'babs1': 170. * qcc.conversion_factor('degree', 'radian'),
                'babs2': 85. * qcc.conversion_factor('degree', 'radian'),
                'babs3': 170. * qcc.conversion_factor('degree', 'radian'),
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
                geo1 = automol.convert.zmatrix.geometry(rct1_zma)
                if not automol.geom.is_linear(geo1):
                    tors_name = ts_name_dct['babs2']
                    tors_names += (tors_name,)

            # babs3 should only be included if there is only group
            # connected to the radical atom
            include = include_babs3(frm_bnd_key, rct2_gra)
            if 'babs3' in ts_name_dct and include:
                tors_name = ts_name_dct['babs3']
                tors_names += (tors_name,)

            frm_bnd_key = shift_vals_from_dummy(frm_bnd_key, ts_zma)
            brk_bnd_key = shift_vals_from_dummy(brk_bnd_key, ts_zma)

            # Build reactants graph
            atom_keys2 = automol.graph.atom_keys(rct2_gra)
            natom2 = len(atom_keys2)
            atm_key_dct = dict(zip(
                atom_keys2,
                (key+natom2 for key in atom_keys2),
            ))
            new_rct2_gra = automol.graph.relabel(rct2_gra, atm_key_dct)
            rcts_gra = automol.graph.union_from_sequence(
                (rct1_gra, new_rct2_gra))

            ret = (ts_zma, dist_name,
                   frm_bnd_key, brk_bnd_key,
                   tors_names, rcts_gra)

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
        rct2_gra = automol.zmatrix.graph(rct_zmas[1], remove_stereo=True)
        rad_atm_keys = automol.graph.resonance_dominant_radical_atom_keys(
            rct2_gra)
        # hno2 hack
        # rad_atm_keys = [0]
        # print('rad_atm_keys:', rad_atm_keys)
        # print('rct_zmas:', rct_zmas)
        # import sys
        # sys.exit()
        if 0 not in rad_atm_keys:
            rct_zmas[1] = reorder_zma_for_radicals(
                rct_zmas[1], min(rad_atm_keys))
            rct2_gra = automol.zmatrix.graph(rct_zmas[1], remove_stereo=True)
            rad_atm_keys = automol.graph.resonance_dominant_radical_atom_keys(
                rct2_gra)
            # following assert checks to ensure that the first atom
            # in the second reactant is a radical
            # this is required for the remainder of the routine
            assert 0 in rad_atm_keys
        rct_zmas, rct_gras = shifted_standard_zmas_graphs(
            rct_zmas, remove_stereo=True)
        prd_zmas, prd_gras = shifted_standard_zmas_graphs(
            prd_zmas, remove_stereo=True)
        # fix to put radical atom first
        # ultimately need to fix this for multiple radical centers
        # end of fix
        # print(rct_gras)
        tras, _, _ = automol.graph.reac.hydrogen_abstraction(
            rct_gras, prd_gras)
        # print('tras')
        # for tra in tras:
        if tras:
            tra = tras[0]
            rct1_gra, rct2_gra = rct_gras
            rct1_zma, rct2_zma = rct_zmas
            rct1_natms = automol.zmatrix.count(rct1_zma)
            rct2_natms = automol.zmatrix.count(rct2_zma)

            frm_bnd_key, = automol.graph.trans.formed_bond_keys(tra)
            brk_bnd_key, = automol.graph.trans.broken_bond_keys(tra)

            rct2_atm1_key = _xor(frm_bnd_key, brk_bnd_key)
            if not rct2_atm1_key == rct1_natms:
                rct2_gra = automol.graph.move_idx_to_top(rct2_gra, rct2_atm1_key, rct1_natms)
                rct2_zma = automol.geom.zmatrix(automol.graph.geometry(rct2_gra))
                rct_zmas = [rct1_zma, rct2_zma]
                rct_zmas, _ = shifted_standard_zmas_graphs(rct_zmas, remove_stereo=True)
                _, rct2_zma = rct_zmas
                new_bnd_key = []    
                for key in frm_bnd_key:
                    if key == rct2_atm1_key:
                        key = rct1_natms
                    new_bnd_key.append(key)
                frm_bnd_key = frozenset(new_bnd_key)    
                rct_gras = (rct1_gra, rct2_gra)
            #print('frm_key', frm_bnd_key)
            #print('brk_key', brk_bnd_key)
            rct1_atm1_key = next(iter(frm_bnd_key & brk_bnd_key))

            # if rct1 and rct2 are isomorphic, we may get an atom key on rct2.
            # in that case, determine the equivalent atom from rct1
            if rct1_atm1_key in automol.graph.atom_keys(rct2_gra):
                atm_key_dct = automol.graph.full_isomorphism(rct2_gra,
                                                             rct1_gra)
                assert atm_key_dct
                rct1_atm1_key = atm_key_dct[rct1_atm1_key]

            rct1_atm2_key, rct1_atm3_key, _ = join_atom_keys(
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
                'babs1': 170. * qcc.conversion_factor('degree', 'radian'),
                'babs2': 85. * qcc.conversion_factor('degree', 'radian'),
                'babs3': 85. * qcc.conversion_factor('degree', 'radian'),
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
                geo1 = automol.convert.zmatrix.geometry(rct1_zma)
                if not automol.geom.is_linear(geo1):
                    tors_name = ts_name_dct['babs2']
                    tors_names += (tors_name,)

            include = include_babs3(frm_bnd_key, rct2_gra)
            if 'babs3' in ts_name_dct and include:
                tors_name = ts_name_dct['babs3']
                tors_names += (tors_name,)

            frm_bnd_key = shift_vals_from_dummy(frm_bnd_key, ts_zma)
            brk_bnd_key = shift_vals_from_dummy(brk_bnd_key, ts_zma)

            # Build reactants graph
            atom_keys2 = automol.graph.atom_keys(rct2_gra)
            natom2 = len(atom_keys2)
            atm_key_dct = dict(zip(
                atom_keys2,
                (key+natom2 for key in atom_keys2),
            ))
            new_rct2_gra = automol.graph.relabel(rct2_gra, atm_key_dct)
            rcts_gra = automol.graph.union_from_sequence(
                (rct1_gra, new_rct2_gra))
    
            # print('frm', frm_bnd_key)
            # print('brk', brk_bnd_key)

            ret = (ts_zma, dist_name,
                   frm_bnd_key, brk_bnd_key,
                   tors_names, rcts_gra)

    return ret

def _xor(atms1, atms2):
    for atmi in atms1:
        if not atmi in atms2:
            return atmi

def shift_vals_from_dummy(vals, zma):
    """ Shift a set of values using remdummy
        Shift requires indices be 1-indexed
    """
    type_ = type(vals)

    dummy_idxs = automol.zmatrix.atom_indices(zma, sym='X')

    shift_vals = []
    for val in vals:
        shift = 0
        for dummy in dummy_idxs:
            if val >= dummy:
                shift += 1
        shift_vals.append(val+shift)

    shift_vals = type_(shift_vals)

    return shift_vals


if __name__ == '__main__':
    import automol

    RCT_ICHS = (
        automol.smiles.inchi('[O]O'),
        automol.smiles.inchi('[CH2]C=CCCCCCCCCC'),
    )
    #PRD_ICHS1 = (
    #    automol.smiles.inchi('[O]O'),
    #    automol.smiles.inchi('[CH2]C=CCC'),
    #)
    #PRD_ICHS2 = (
    #    automol.smiles.inchi('[O]O'),
    #    automol.smiles.inchi('CC=C[CH]C'),
    #)
    PRD_ICHS1 = (
        automol.smiles.inchi('[O][O]'),
        automol.smiles.inchi('C=CCCCCCCCCCC'),
    )
    PRD_ICHS2 = (
        automol.smiles.inchi('[O][O]'),
        automol.smiles.inchi('CC=CCCCCCCCCC'),
    )
    RCT_ZMAS = [automol.geom.zmatrix(automol.inchi.geometry(ich)) for ich in RCT_ICHS]
    PRD1_ZMAS = [automol.geom.zmatrix(automol.inchi.geometry(ich)) for ich in PRD_ICHS1]
    PRD2_ZMAS = [automol.geom.zmatrix(automol.inchi.geometry(ich)) for ich in PRD_ICHS2]

    RCT_GRAS, _ = automol.graph.standard_keys_for_sequence([
        automol.graph.explicit(automol.inchi.graph(ich, no_stereo=True))
        for ich in RCT_ICHS
    ])
    PRD_GRAS1, _ = automol.graph.standard_keys_for_sequence([
        automol.graph.explicit(automol.inchi.graph(ich, no_stereo=True))
        for ich in PRD_ICHS1
    ])
    PRD_GRAS2, _ = automol.graph.standard_keys_for_sequence([
        automol.graph.explicit(automol.inchi.graph(ich, no_stereo=True))
        for ich in PRD_ICHS2
    ])

    #TS1 = hydrogen_abstraction(RCT_GRAS, PRD_GRAS1)
    #TS2 = hydrogen_abstraction(RCT_GRAS, PRD_GRAS2)
    #TS1 = hydrogen_abstraction(PRD_GRAS1, RCT_GRAS)
    #TS2 = hydrogen_abstraction(PRD_GRAS2, RCT_GRAS)
    TS1 = hydrogen_abstraction(RCT_ZMAS, PRD1_ZMAS)
    TS2 = hydrogen_abstraction(RCT_ZMAS, PRD2_ZMAS)
    #print(automol.geom.string(automol.zmatrix.geometry(TS1[0])))
    #print()
    #print(automol.geom.string(automol.zmatrix.geometry(TS2[0])))
