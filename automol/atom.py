""" functions operating on atomic symbols
"""

SYMBOLS = ['X',
           'H', 'He',
           'C',
           'N',
           'O', 'S',
           'F', 'Cl',
           'Ne', 'Ar']


def standard_case(sym):
    """ returns a camel-case atomic symbol
    """
    return str(sym).title()


def valence(sym):
    """ bonding valence
    """
    return {'X': 0,
            'H': 1, 'He': 0,
            'C': 4,
            'N': 3,
            'O': 2, 'S': 2,
            'F': 1, 'Cl': 1,
            'Ne': 0, 'Ar': 0}[standard_case(sym)]


def nuclear_charge(sym):
    """ nuclear charge
    """
    return {'X': 0,
            'H': 1, 'He': 2,
            'C': 6,
            'N': 7,
            'O': 8, 'S': 16,
            'F': 9, 'Cl': 17,
            'Ne': 10, 'Ar': 18}[standard_case(sym)]


# def lone_pair_count(sym):
#     """ lone pair count
#     """
#     return {'X': 0,
#             'H': 0, 'HE': 1,
#             'C': 0,
#             'N': 1,
#             'O': 2, 'S': 2,
#             'F': 3, 'CL': 3,
#             'NE': 4, 'AR': 4}[sym.upper()]
