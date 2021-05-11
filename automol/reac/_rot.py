""" rotational bond/torsion info for specific reaction classes
"""
from automol.par import ReactionClass
import automol.zmat
from automol.reac._util import hydrogen_abstraction_atom_keys
from automol.reac._util import substitution_atom_keys


# Bimolecular reactions
# 1. Hydrogen abstractions
def hydrogen_abstraction_linear_atom_keys(rxn, zma=None):
    """ Obtain the linear atom keys for a hydrogen abstraction

    :param rxn: a Reaction object
    :param zma: a z-matrix; if passed in, the linear atoms will be determined
        from this; otherwise they will be determined heuristically from the
        reaction object
    :returns: the keys of the linear atoms in the graph
    :rtype: tuple[int]
    """
    tsg = rxn.forward_ts_graph
    if zma is not None:
        lin_keys = list(automol.zmat.linear_atom_keys(zma))
    else:
        lin_keys = list(automol.graph.linear_atom_keys(tsg))

    _, hyd_key, _ = hydrogen_abstraction_atom_keys(rxn)

    lin_keys.append(hyd_key)
    lin_keys = tuple(sorted(set(lin_keys)))
    return lin_keys


# 4. Substitution
def substitution_linear_atom_keys(rxn, zma=None):
    """ Obtain the linear atom keys for a substitution

    :param rxn: a Reaction object
    :param zma: a z-matrix; if passed in, the linear atoms will be determined
        from this; otherwise they will be determined heuristically from the
        reaction object
    :returns: the keys of the linear atoms in the graph
    :rtype: tuple[int]
    """
    tsg = rxn.forward_ts_graph
    if zma is not None:
        lin_keys = list(automol.zmat.linear_atom_keys(zma))
    else:
        lin_keys = list(automol.graph.linear_atom_keys(tsg))

    _, tra_key, _ = substitution_atom_keys(rxn)

    lin_keys.append(tra_key)
    lin_keys = tuple(sorted(set(lin_keys)))
    return lin_keys


def linear_atom_keys(rxn, zma=None):
    """ Obtain the linear atom keys

    :param rxn: a hydrogen migration Reaction object
    :param zma: a z-matrix; if passed in, the linear atoms will be determined
        from this; otherwise they will be determined heuristically from the
        reaction object
    :returns: the keys of the linear atoms in the graph
    :rtype: tuple[int]
    """

    def _default(rxn, zma=None):
        tsg = rxn.forward_ts_graph
        if zma is not None:
            lin_keys = automol.zmat.linear_atom_keys(zma)
        else:
            lin_keys = automol.graph.linear_atom_keys(tsg)
        return lin_keys

    function_dct = {
        # unimolecular
        ReactionClass.HYDROGEN_MIGRATION: _default,
        ReactionClass.BETA_SCISSION: _default,
        ReactionClass.RING_FORM_SCISSION: _default,
        ReactionClass.ELIMINATION: _default,
        # bimolecular
        ReactionClass.HYDROGEN_ABSTRACTION:
        hydrogen_abstraction_linear_atom_keys,
        ReactionClass.ADDITION: _default,
        ReactionClass.INSERTION: _default,
        ReactionClass.SUBSTITUTION: substitution_linear_atom_keys,
    }

    fun_ = function_dct[rxn.class_]
    ret = fun_(rxn, zma=zma)
    return ret


def rotational_bond_keys(rxn, zma=None):
    """ Obtain the rotational bond keys

    :param rxn: a hydrogen migration Reaction object
    :param zma: a z-matrix; if passed in, the linear atoms will be determined
        from this; otherwise they will be determined heuristically from the
        reaction object
    :returns: the keys of the rotational bonds in the graph
    :rtype: tuple[frozenset[int]]
    """
    tsg = rxn.forward_ts_graph
    lin_keys = linear_atom_keys(rxn, zma=zma)
    bnd_keys = automol.graph.rotational_bond_keys(tsg, lin_keys=lin_keys)
    return bnd_keys


def rotational_groups(rxn, key1, key2, dummy=False):
    """ Obtain the rotational groups for a given rotational axis

    :param rxn: a hydrogen migration Reaction object
    :param zma: a z-matrix; if passed in, the linear atoms will be determined
        from this; otherwise they will be determined heuristically from the
        reaction object
    :returns: the rotational groups on either side of the axis
    :rtype: (tuple[int], tuple[int])
    """
    tsg = rxn.forward_ts_graph
    grps = automol.graph.rotational_groups(tsg, key1, key2, dummy=dummy)
    return grps


def rotational_symmetry_number(rxn, key1, key2, zma=None):
    """ Obtain the rotational symmetry number for a given rotational axis

    :param rxn: a hydrogen migration Reaction object
    :param zma: a z-matrix; if passed in, the linear atoms will be determined
        from this; otherwise they will be determined heuristically from the
        reaction object
    :returns: the rotational symmetry number of the axis
    :rtype: int
    """
    lin_keys = linear_atom_keys(rxn, zma=zma)
    tsg = rxn.forward_ts_graph
    sym_num = automol.graph.rotational_symmetry_number(tsg, key1, key2,
                                                       lin_keys=lin_keys)
    return sym_num
