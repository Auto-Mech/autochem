"""
  Various info objects for torsions

  maybe build rotor dcts here since they are just collections of torsions?
"""

import numpy
import yaml
import automol
from phydat import phycon


class Parameter():
    """ Info which can define a new torsion (need new class name)
    """

    ROT_GROUP = 'group'
    SYMMETRY = 'symmetry',
    POTENTIAL = 'pot'
    SPAN = 'span'
    GEOMETRY = 'geo'
    ZMATRIX = 'zma'


# constructor
def from_zmatrix(zma, tors_name,
                 scan_increment=30.0, tors_model='1dhr',
                 frm_bnd_keys=(), brk_bnd_keys=()):
    """ Build the rotor object consisting of all of the torsions
        and their info.

        :param zma: zmatrix
        :type zma:
        :param rotor_names: names of torsions in rotors
            [ ['D1'], ['D2', 'D3'], ['D4'] ]
        :type rotor_names: tuple(tuple(str))
    """

    tors_dct = {
        Parameter.ROT_GROUP: torsional_groups(),
        Parameter.SYMMETRY; symmetry_number(),
        Parameter.SPAN: span(),
        Parameter.GRID: grid()
        Parameter.GRID: potential()
    }

    # Check symmetry number
    grid = torsional_grids()

    # Structures of torsion
    # mode_idxs = automol.zmatrix.coord_idxs(zma, tname)
    # mode_idxs = tuple((idx+1 for idx in mode_idxs))

    return rotor_dct


# getters
def symmetry_number(tors):
    """ Get the symmery number of a torsion.

        :param tors: torsion object
        :type tors: dict[str: obj]
        :rtype: int
    """
    return tors.get(Parameter.SYMMETRY, None)


def grid(tors):
    """ Get the grid of a torsion.

        :param tors: torsion object
        :type tors: dict[str: obj]
        :rtype: int
    """
    return tors.get(Parameter.GRID, None)


def span(tors):
    """ Get the span of a torsion.

        :param tors: torsion object
        :type tors: dict[str: obj]
        :rtype: float
    """
    return tors.get(Parameter.SPAN, None)


def axis(tors):
    """ Get the axis of a torsion.

        :param tors: torsion object
        :type tors: dict[str: obj]
        :rtype: (int, int)
    """

    axis1 = tors.get(Parameter.AXIS1, None)
    axis2 = tors.get(Parameter.AXIS2, None)
    if axis1 is not None and axis2 is not None:
        _axis = (axis1, axis2)
    else:
        _axis = None

    return _axis


def rotational_groups(tors):
    """ Get the rotational groups of a torsion.

        :param tors: torsion object
        :type tors: dict[str: obj]
        :rtype: (tuple(int), tuple(int))
    """

    group1 = tors.get(Parameter.GROUP1, None)
    group2 = tors.get(Parameter.GROUP2, None)
    if group1 is not None and group2 is not None:
        _group = (group1, group2)
    else:
        _group = None

    return _group


def geometry(tors):
    """ Get the molecular geometry associated with torsion.

        :param tors: torsion object
        :type tors: dict[str: obj]
        :rtype: int
    """
    return tors.get(Parameter.GEO, None)


def zmatrix(tors):
    """ Get the Z-Matrix associated with torsion.

        :param tors: torsion object
        :type tors: dict[str: obj]
        :rtype: int
    """
    return tors.get(Parameter.ZMA, None)
