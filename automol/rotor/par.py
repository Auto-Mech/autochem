"""
    Torsional parameters
"""

class Torsion():
    """ Info which can define a new torsion (need new class name)
    """

    GROUP1 = 'group1'
    GROUP2 = 'group2'
    AXIS1 = 'axis1'
    AXIS2 = 'axis2'
    SYMMETRY = 'symmetry',
    POTENTIAL = 'pot'
    SPAN = 'span'
    GEOMETRY = 'geo'
    ZMATRIX = 'zma'


class Model():
    """ Torsional models
    """
    1DHR = '1dhr'
    1DHR_FRZ_TORS = '1dhrf'
    1DHR_FRZ_ALL = '1dhrfa'
    MDHR = 'mdhr'
    MDHR_VIB_ADIAB = 'mdhrv'
