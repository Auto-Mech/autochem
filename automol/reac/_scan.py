""" scan info for specific reaction classes

(Currently just the scan coordinate, but we may want to add in other things)
"""
import more_itertools as mit
from automol.graph import ts
from automol import par
import automol.zmat
from automol.reac._util import hydrogen_migration_atom_keys
from automol.reac._util import ring_forming_scission_chain
from automol.reac._util import insertion_forming_bond_keys


# Unimolecular reactions
# 1. Hydrogen migrations
def hydrogen_migration_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for a hydrogen migration.

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return scan_name


def hydrogen_migration_constraint_coordinates(rxn, zma):
    """ Obtain the constraint coordinates for a hydrogen migration

    :param rxn: a Reaction object
    :returns: the names of the constraint coordinates in the z-matrix
    :rtype: tuple[str]
    """
    att_key, _, _, ngb_key = hydrogen_migration_atom_keys(rxn)
    dist_name = automol.zmat.distance_coordinate_name(zma, att_key, ngb_key)
    return (dist_name,)


# 2. Beta scissions
def beta_scission_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for a beta scission

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *brk_bnd_key)
    return scan_name


# 3. Ring-forming scissions
def ring_forming_scission_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for a ring-forming scission

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *brk_bnd_key)
    return scan_name


def ring_forming_scission_constraint_coordinates(rxn, zma):
    """ Obtain the constraint coordinates for a ring-forming scission

    :param rxn: a Reaction object
    :returns: the names of the constraint coordinates in the z-matrix
    :rtype: str
    """
    chain_keys = ring_forming_scission_chain(rxn)
    ang_keys_lst = sorted(mit.windowed(chain_keys[1:], 3))
    dih_keys_lst = sorted(mit.windowed(chain_keys, 4))
    ang_names = [automol.zmat.central_angle_coordinate_name(zma, *ks)
                 for ks in ang_keys_lst]
    dih_names = [automol.zmat.dihedral_angle_coordinate_name(zma, *ks)
                 for ks in dih_keys_lst]
    const_names = tuple(ang_names + dih_names)
    return const_names


# 4. Eliminations
def elimination_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for an elimination

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return scan_name


# Bimolecular reactions
# 1. Hydrogen abstractions
def hydrogen_abstraction_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for a hydrogen abstraction

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return scan_name


# 2. Additions
def addition_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for an addition

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return scan_name


# 3. Insertions
def insertion_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for an insertion

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, _ = insertion_forming_bond_keys(rxn)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return scan_name


# 4. Substitution
def substitution_scan_coordinate(rxn, zma):
    """ Obtain the scan coordinate for a substitution

    :param rxn: a Reaction object
    :returns: the name of the scan coordinate in the z-matrix
    :rtype: str
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    scan_name = automol.zmat.distance_coordinate_name(zma, *frm_bnd_key)
    return scan_name


def scan_coordinate(rxn, zma):
    """ Obtain the scan coordinates

    :param rxn: a hydrogen migration Reaction object
    """
    function_dct = {
        # unimolecular
        par.ReactionClass.HYDROGEN_MIGRATION:
        hydrogen_migration_scan_coordinate,
        par.ReactionClass.BETA_SCISSION: beta_scission_scan_coordinate,
        par.ReactionClass.RING_FORM_SCISSION:
        ring_forming_scission_scan_coordinate,
        par.ReactionClass.ELIMINATION: elimination_scan_coordinate,
        # bimolecular
        par.ReactionClass.HYDROGEN_ABSTRACTION:
        hydrogen_abstraction_scan_coordinate,
        par.ReactionClass.ADDITION: addition_scan_coordinate,
        par.ReactionClass.INSERTION: insertion_scan_coordinate,
        par.ReactionClass.SUBSTITUTION: substitution_scan_coordinate,
    }

    fun_ = function_dct[rxn.class_]
    ret = fun_(rxn, zma)
    return ret


def constraint_coordinates(rxn, zma):
    """ Obtain the constraint coordinates

    :param rxn: a hydrogen migration Reaction object
    """

    def _return_empty_tuple(*_):
        return ()

    function_dct = {
        # unimolecular
        par.ReactionClass.HYDROGEN_MIGRATION:
        hydrogen_migration_constraint_coordinates,
        par.ReactionClass.BETA_SCISSION: _return_empty_tuple,
        par.ReactionClass.RING_FORM_SCISSION:
        ring_forming_scission_constraint_coordinates,
        par.ReactionClass.ELIMINATION: _return_empty_tuple,
        # bimolecular
        par.ReactionClass.HYDROGEN_ABSTRACTION:
        _return_empty_tuple,
        par.ReactionClass.ADDITION: _return_empty_tuple,
        par.ReactionClass.INSERTION: _return_empty_tuple,
        par.ReactionClass.SUBSTITUTION: _return_empty_tuple,
    }

    fun_ = function_dct[rxn.class_]
    ret = fun_(rxn, zma)
    return ret



