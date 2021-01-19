""" scan info for specific reaction classes

(Currently just the scan coordinate, but we may want to add in other things)
"""
from automol.graph import ts
from automol import par


# Unimolecular reactions
# 1. Hydrogen migrations
def hydrogen_migration_scan_coordinate(rxn):
    """ scan coordinate for a hydrogen migration

    :param rxn: a Reaction object
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    return frm_bnd_key


# 2. Beta scissions
def beta_scission_scan_coordinate(rxn):
    """ scan coordinate for a beta scission

    :param rxn: a Reaction object
    """
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    return brk_bnd_key


# 3. Ring-forming scissions
def ring_forming_scission_scan_coordinate(rxn):
    """ scan coordinate for a ring-forming scission

    :param rxn: a Reaction object
    """
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    return brk_bnd_key


# 4. Eliminations
def elimination_scan_coordinate(rxn):
    """ scan coordinate for an elimination

    :param rxn: a Reaction object
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    return frm_bnd_key


# Bimolecular reactions
# 1. Hydrogen abstractions
def hydrogen_abstraction_scan_coordinate(rxn):
    """ scan coordinate for a hydrogen abstraction

    :param rxn: a Reaction object
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    return frm_bnd_key


# 2. Additions
def addition_scan_coordinate(rxn):
    """ scan coordinate for an addition

    :param rxn: a Reaction object
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    return frm_bnd_key


# 3. Insertions
def insertion_scan_coordinate(rxn):
    """ scan coordinate for an insertion

    :param rxn: a Reaction object
    """
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    # Choose the forming bond that doesn't intersect with the breaking bond, if
    # one of them does
    frm_bnd_keys = sorted(ts.forming_bond_keys(rxn.forward_ts_graph),
                          key=sorted)
    frm_bnd_keys = sorted(frm_bnd_keys,
                          key=lambda x: len(x & brk_bnd_key))
    frm_bnd_key = frm_bnd_keys[0]
    return frm_bnd_key


# 4. Substitution
def substitution_scan_coordinate(rxn):
    """ scan coordinate for a substitution

    :param rxn: a Reaction object
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    return frm_bnd_key


def scan_coordinate(rxn):
    """ scan coordinates

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
    ret = fun_(rxn)
    return ret
