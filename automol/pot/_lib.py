"""
  Library of parameters
"""


# DCTS OF POTENTIAL PARAMETERS
# eps[whatever], sig[ang] params
LJ_DCT = {
    ('H', 'H'): (0.25, 1.0),
    ('H', 'C'): (0.25, 1.0),
    ('H', 'O'): (0.25, 1.0),
    ('C', 'C'): (0.25, 1.0),
    ('C', 'O'): (0.25, 1.0),
    ('O', 'O'): (0.25, 1.0)
}

# A, B, C params E[kcal] R[Ang]; R cutoff
EXP6_DCT = {
    ('H', 'H'): (2.442e3, 3.74, 48.8, 1.0),
    ('H', 'C'): (6.45e3, 3.67, 116.0, 1.0),
    ('H', 'N'): (6.45e3, 3.67, 116.0, 1.0),
    ('H', 'O'): (6.45e3, 3.67, 116.0, 1.0),
    ('C', 'C'): (7.69e4, 3.6, 460.0, 0.8),
    ('C', 'O'): (7.69e4, 3.6, 460.0, 0.8),
    ('O', 'O'): (7.69e4, 3.6, 460.0, 0.8),
    ('Cl', 'Cl'): (0.0, 3.6, 460.0, 0.8),
    ('Cl', 'C'): (0.0, 3.6, 460.0, 0.8),
    ('Cl', 'O'): (0.0, 3.6, 460.0, 0.8),
    ('Cl', 'H'): (0.0, 3.6, 460.0, 0.8),
    ('N', 'N'): (7.69e4, 3.6, 460.0, 0.8),
    ('N', 'O'): (7.69e4, 3.6, 460.0, 0.8),
    ('N', 'C'): (7.69e4, 3.6, 460.0, 0.8)
}


def read_params(dct, symb1, symb2):
    """ Read the parameters from one of the potential
        parameter dictionaries

        :param dct: potential parameter dct
        :type dct: dict[(symbs):(params)]
        :param symb1: atomic symbol of atom1 involved in the interaction
        :param symb1: str
        :param symb2: atomic symbol of atom2 involved in the interaction
        :param symb2: str
        :rtype: tuple(float)
    """

    params = dct.get((symb1, symb2), None)
    if params is None:
        params = dct.get((symb2, symb1), None)

    return params
