"""
  Various info objects for torsions
"""

import numpy
import yaml
import automol
from phydat import phycon


# z-matrix torsional degrees of freedom
def symmetry_number(zma, tors_name,
                    frm_bnd_keys=None, brk_bnd_keys=None):
    """ symmetry numbers for torsional dihedrals
    """

    dih_edg_key_dct = _dihedral_edge_keys(zma)
    assert {tors_name} <= set(dih_edg_key_dct.keys())
    edg_key = dih_edg_key_dct[tors_name]

    gra = automol.convert.zmatrix.graph(zma, remove_stereo=True)
    bnd_sym_num_dct = automol.graph.bond_symmetry_numbers(
        gra, frm_bnd_keys, brk_bnd_keys)

    if edg_key in bnd_sym_num_dct.keys():
        sym_num = bnd_sym_num_dct[edg_key]
    else:
        sym_num = 1

    return sym_num


def grid(zma, tors_name, scan_increment, frm_bnd_keys=None, brk_bnd_keys=None):
    """ scan grids
    """

    val_dct = automol.zmat.value_dictionary(zma)

    linspace = _scan_linspace(
        zma, tors_name, scan_increment=scan_increment,
        frm_bnd_keys=frm_bnd_keys, brk_bnd_keys=brk_bnd_keys)
    _grid = numpy.linspace(*linspace)

    ini_val = val_dct[tors_name]
    grid_from_equil = tuple(val.item() + ini_val for val in _grid)

    return grid_from_equil


def span(zma, tors_name, frm_bnd_keys=None, brk_bnd_keys=None):
    """ Determine the span of the torsion
    """
    sym_num = symmetry_number(
        zma, tors_name,
        frm_bnd_keys=frm_bnd_keys, brk_bnd_keys=brk_bnd_keys)
    return (2.0 * numpy.pi) / sym_num


# Rotor indices
# Build the remdummy list here and shift things
def axis(zma, tors_name):
    """ Determine the axis about which the rotational groups
        of the torsion group rotate: R1-X1-X2-R2
    """
    coo_dct = automol.zmat.coordinates(zma, multi=False)
    return coo_dct[tors_name][1:3]


def torsional_groups(zma, tors_name):
    """ Build the rotational groups (need to add TS mode)
    """

    # Build the axis
    tors_axis = axis(zma, tors_name)
    tors_key, ts_bnd = _tors_key(tors_axis, frm_bnd_key=None, brk_bnd_key=None)

    # Initialize the group
    gra = automol.zmat.graph(zma, remove_stereo=True)
    if ts_bnd not in automol.graph.bond_keys(gra):
        gra = automol.graph.add_ts_bonds(gra, keys=[ts_bnd])
    group = list(
        automol.graph.branch_atom_keys(
            gra, tors_key, axis) - set(axis))
    if not group:
        for atm in tors_axis:
            if atm != tors_key:
                tors_key = atm
        group = list(
            automol.graph.branch_atom_keys(
                gra, tors_key, axis) - set(axis))

    # augment the groups as needed
    _add_ts_bnd_to_torsional_group(zma, ts_bnd, rxn_class=None)

    # next augment


def _tors_key(tors_axis, frm_bnd_keys=None, brk_bnd_keys=None):
    """ build a key for a bond
    """

    # Set key to initial guess
    tors_key = tors_axis[1]

    # Re-set key for a transition state
    if frm_bnd_keys:
        ts_bnd = frm_bnd_keys
    elif brk_bnd_keys:
        ts_bnd = brk_bnd_keys
    else:
        ts_bnd = ()

    if ts_bnd:
        for atm in tors_axis:
            if atm in ts_bnd:
                tors_key = atm
                break

    return tors_key, ts_bnd


def _add_ts_bnd_to_torsional_group(zma, ts_bnd, rxn_class=None):
    """ Add the group for a torsion to rotational group
    """
    n_atm = automol.zmat.count(zma)
    if 'addition' in rxn_class or 'abstraction' in rxn_class:
        group2 = []
        ts_bnd1 = min(ts_bnd)
        ts_bnd2 = max(ts_bnd)
        for idx in range(ts_bnd2, n_atm):
            group2.append(idx)
        if ts_bnd1 in group:
            for atm in group2:
                if atm not in group:
                    group.append(atm)

    return group

# remdummy = geomprep.build_remdummy_shift_lst(zma)

# I/O
def string(tors_dct):
    """ write the transformation to a string
    """

    def _encode_group(group_keys):
        group_str = '-'.join(str(val) for val in group_keys)
        return group_str

    tors_names = _sort_tors_names(list(tors_dct.keys()))

    str_dct = {}
    for name in tors_names:
        dct = tors_dct[name]
        str_dct[name] = {
            'axis': _encode_group(dct.get('axis', None)),
            'group1': _encode_group(dct.get('group1', None)),
            'group2': _encode_group(dct.get('group2', None)),
            'symmetry': dct.get('symmetry', None),
            'span': round(dct.get('span', None) * phycon.RAD2DEG, 2)
        }

    tors_str = yaml.dump(str_dct, sort_keys=False)

    return tors_str


def from_string(tors_str):
    """ read the transformation from a string
    """

    def _decode_group(group_str):
        group_keys = tuple(map(int, group_str.split('-')))
        return group_keys

    tors_dct = yaml.load(tors_str, Loader=yaml.FullLoader)
    for dct in tors_dct.values():
        dct['axis'] = _decode_group(dct['axis'])
        dct['group1'] = _decode_group(dct['group1'])
        dct['group2'] = _decode_group(dct['group2'])
        dct['span'] = dct['span'] * phycon.DEG2RAD

    return tors_dct


# Helper functions
def _dihedral_edge_keys(zma):
    """ dihedral bonds, by name
    """

    coo_dct = automol.zmat.coordinates(zma)
    dih_names = automol.zmat.dihedral_angle_names(zma)
    dih_keys_lst = tuple(map(coo_dct.__getitem__, dih_names))
    dih_edg_key_dct = {dih_name: frozenset(dih_key[1:3])
                       for dih_name, dih_keys in zip(dih_names, dih_keys_lst)
                       for dih_key in dih_keys}

    return dih_edg_key_dct


def _scan_linspace(zma, tors_name, scan_increment=0.523599,
                   frm_bnd_keys=None, brk_bnd_keys=None):
    """ scan grids for torsional dihedrals
    """

    sym_num = symmetry_number(
        zma, tors_name,
        frm_bnd_keys=frm_bnd_keys, brk_bnd_keys=brk_bnd_keys)
    interval = ((2.0 * numpy.pi) / sym_num) - scan_increment
    npoints = int(interval / scan_increment) + 2

    return (0.0, interval, npoints)


def _sort_tors_names(tors_names):
    """ sort torsional names
    """
    tors_names = list(tors_names)
    tors_names.sort(key=lambda x: int(x.split('D')[1]))
    return tors_names
