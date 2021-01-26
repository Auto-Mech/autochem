"""
  Functions handling hindered rotor model calculations

    rotor_dct = {
        (tors_name1, tors_name2, tors_name3): ({tors_dct1}, {tors_dct2}, {tors_dct3}), {inf_dct #pot})
        # put a zma/geo with it


    rotors = (
        rotor_dct1 = {
            rotor_tors_name1: {tors_inf: tors_bj},
            rotor_tors_name2: {tors_inf: tors_bj},
        },
        rotor_dct2 = {
            rotor_tors_name1: {tors_inf: tors_bj},
            rotor_tors_name2: {tors_inf: tors_bj},
        }
    )

"""


def from_tors_dcts(tors_dcts):
    """ Build 
    """


def from_zma(zma):
    """
    """

    # Determine all of the torsional names
    tors_names = []

    # Build all of the torsional dcts
    tors_dcts = {}

    # Combine into list of 1DHR rotors
    rotors = (
        {tors_name1: tors_dct1},
        {tors_name2: tors_dct2}
    )

    # Recombine into multirotors if possible
    # probably best to just combine everything into one multirotor and then 

    # Group the torsional names into rotor combinations
    rotor_names = ()

    for name in tors_names


def rebuild(rotor
