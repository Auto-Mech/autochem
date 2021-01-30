""" Handle name groups
"""

import automol


def group_torsions_into_rotors(tors_lst, name_grps=None, tors_model=None):
    """ take a list of torsion objects and build the rotors

        Uses the names to build the the lists
    """

    tors_names = tuple(tors.name for tors in tors_lst)
    name_tors_dct = dict(zip(tors_names, tors_lst))

    if name_grps is None:
        name_grps = name_groups(tors_names, model='1dhr')

    rotors = ()
    for name_grp in name_grps:
        rotor = ()
        for name in name_grp:
            rotor += (name_tors_dct[name],)
        rotors += (rotor,)

    return rotors


def name_groups(names, model):
    """ Build the list of names
    """

    if model in ('1dhr', '1dhrf', '1dhrfa', 'tau'):
        tors_names = [[name] for name in names]
    else:
        tors_names = [[name for name in names]]

    tors_names = tuple(tuple(x) for x in tors_names)

    return tors_names


def mdhr_prep(zma, run_tors_names):
    """ Handle cases where the MDHR
    """

    # Figure out set of torsions are to be used: defined or AMech generated
    rotor_lst = run_tors_names

    # Check the dimensionality of each rotor to see if they are greater than 4
    # Call a function to reduce large rotors
    final_rotor_lst = []
    for torsion in rotor_lst:
        if len(rotor) > 4:
            for reduced_rotor in reduce_rotor_dimensionality(zma, rotor):
                final_rotor_lst.append(reduced_rotor)
        else:
            final_rotor_lst.append(rotor)

    return final_rotor_lst


def reduce_rotor_dimensionality(zma, rotor):
    """ For rotors with a dimensionality greater than 4, try and take them out
    """

    reduced_rotor_lst = []
    methyl_rotors = []
    for tors in rotor:
        # If a methyl rotor add to methyl rotor list, or add to reduced lst
        if is_methyl_rotor(zma, rotor):   # Add arguments when ID methyls
            methyl_rotors.append(zma, tors)
        else:
            reduced_rotor_lst.append(tors)

    # Add each of methyl rotors, if any exist
    if methyl_rotors:
        for methyl_rotor in methyl_rotors:
            reduced_rotor_lst.append(methyl_rotor)

    # Check new dimensionality of list; if still high, flatten to lst of 1DHRs
    if len(reduced_rotor_lst) > 4:
        reduced_rotor_lst = [tors
                             for rotor in reduced_rotor_lst
                             for tors in rotor]

    return reduced_rotor_lst


def _is_methyl_rotor(zma, axis, grp1, grp2):
    """ Identify if the rotor
    """

    symbs = automol.zmat.symbols(zma)
    grp1_symbs = sorted(list(symbs[idx] for idx in grp1))
    grp2_symbs = sorted(list(symbs[idx] for idx in grp2))

    return bool(
        any(grp == ['C', 'H', 'H', 'H'] for grp in (grp1_symbs, grp2_symbs)))
