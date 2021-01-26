""" handle multidimensional rotors
"""


def mdhr_prep(zma, run_tors_names):
    """ Handle cases where the MDHR
    """

    # Figure out set of torsions are to be used: defined or AMech generated
    rotor_lst = run_tors_names

    # Check the dimensionality of each rotor to see if they are greater than 4
    # Call a function to reduce large rotors
    final_rotor_lst = []
    for rotor in rotor_lst:
        if len(rotor) > 4:
            for reduced_rotor in reduce_rotor_dimensionality(zma, rotor):
                final_rotor_lst.append(reduced_rotor)
        else:
            final_rotor_lst.append(rotor)

    return final_rotor_lst


def reduce_rotor_dimensionality(zma, rotor):
    """ For rotors with a dimensionality greater than 4, try and take them out
    """

    # Find the methyl rotors for that are a part of the MDHR
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
