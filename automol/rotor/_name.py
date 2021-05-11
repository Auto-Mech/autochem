""" Handle name groups
"""

from automol.zmat import symbols


def group_torsions_into_rotors(tors_lst, name_grps=None, multi=False):
    """ take a list of torsion objects and build the rotors

        Uses the names to build the the lists
    """

    tors_names = tuple(tors.name for tors in tors_lst)
    name_tors_dct = dict(zip(tors_names, tors_lst))

    if name_grps is None:
        name_grps = _name_groups(tors_names, multi=multi)

    rotors = ()
    for name_grp in name_grps:
        rotor = ()
        for name in name_grp:
            rotor += (name_tors_dct[name],)
        rotors += (rotor,)

    # Split up the mdhrs if they are too big
    if multi:
        rotors = _assess_dimensionality(rotors)

    return rotors


def _name_groups(names, multi=False):
    """ Build the list of names
    """

    if not multi:
        tors_names = tuple((name,) for name in names)
    else:
        tors_names = (tuple(name for name in names),)

    return tors_names


def _assess_dimensionality(rotor_lst):
    """ Handle cases where the MDHR
    """

    # Check the dimensionality of each rotor to see if they are greater than 4
    final_rotor = ()
    for rotor in rotor_lst:
        if len(rotor) > 4:
            final_rotor += _reduce_rotor_dimensionality(rotor)
        else:
            final_rotor += (rotor,)

    return final_rotor


def _reduce_rotor_dimensionality(rotor):
    """ For rotors with a dimensionality greater than 4, try and take them out
    """

    # Get the indices of all the -CH3 rotors and non-CH3 rotors
    methyl_idxs = ()
    for idx, tors in enumerate(rotor):
        if _is_methyl_rotor(tors):
            methyl_idxs += (idx,)

    non_methyl_idxs = tuple(i for i in range(len(rotor))
                            if i not in methyl_idxs)

    # Assess if ndim(reduce rotor) > 4; if yes flatten rotor instead of rebuild
    if len(non_methyl_idxs) > 4:
        reduced_rotor = ()
        for tors in rotor:
            reduced_rotor += ((tors,),)
    else:
        reduced_rotor = (tuple(rotor[idx] for idx in non_methyl_idxs),)
        for idx in methyl_idxs:
            reduced_rotor += ((rotor[idx],),)

    return reduced_rotor


def _is_methyl_rotor(tors):
    """ Identify if the rotor
    """

    symbs = symbols(tors.zma)
    axis = list(tors.axis)
    grps = tors.groups
    rgrp1 = [axis[0]] + list(grps[0])
    rgrp2 = [axis[1]] + list(grps[1])
    rgrp1_symbs = sorted(list(symbs[idx] for idx in rgrp1))
    rgrp2_symbs = sorted(list(symbs[idx] for idx in rgrp2))

    return bool(
        any(grp == ['C', 'H', 'H', 'H'] for grp in (rgrp1_symbs, rgrp2_symbs))
    )
