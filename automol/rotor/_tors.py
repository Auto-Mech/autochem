"""
 Build the rotor and torsion objects

 Rotor: (tors_obj_1, tors_obj_2, tors_obj_3)
"""

# import numpy
import yaml
import automol
# from automol.rotor._par import TorsionParam
from phydat import phycon


class Parameter():
    """ Info which can define a new torsion (need new class name)
    """
    GROUPS = 'groups'
    AXIS = 'axis'
    SYMMETRY = 'symmetry'
    POTENTIAL = 'pot'
    SPAN = 'span'
    GEOM = 'geo'
    ZMAT = 'zma'


# class Model():
#     """ Torsional models
#     """
#     ONEDHR = '1dhr'
#     ONEDHR_FRZ_TORS = '1dhrf'
#     1DHR_FRZ_ALL = '1dhrfa'
#     MDHR = 'mdhr'
#     MDHR_VIB_ADIAB = 'mdhrv'


class Torsion:
    """ Describes a torsion, which one or more make up a rotor
    """

    def __init__(self, zma, name, axis, grps, symm):
        """ constructor
        """

        self.name = name
        self.zma = zma
        self.symmetry = symm
        self.groups = grps
        self.axis = axis
        # self.span = span(symmetry)
        # self.indices = Torsion._indices(zma, name)

        # Attributes defaulted to none
        self.pot = None
        self.grid = None

    # @staticmethod
    # def span(symmetry):
    #     """ Obtain the torsional span
    #     """
    #     return (2.0 * numpy.pi) / symmetry

    # @staticmethod
    # def _indices(zma, name):
    #     """ Build indices for the torsion
    #     """
    #     mode_idxs = automol.zmat.coord_idxs(zma, name)
    #     mode_idxs = tuple((idx+1 for idx in mode_idxs))


# Build a converter from object to a dictionary


def torsion_dct(zma, gra, lin_keys):
    """  Build a dictionary of torsion objects indexed by name
    """

    # Build the torsion objects
    _name_axis_dct = name_axis_dct(zma, gra, lin_keys)

    tors_obj_dct = {}
    for name, axis in _name_axis_dct.items():
        grps = groups(gra, axis)
        symm = symmetry(gra, axis, lin_keys)
        tors_obj_dct[name] = Torsion(zma, name, axis, grps, symm)

    return tors_obj_dct


def name_axis_dct(zma, gra, lin_keys):
    """ Generate the bond keys for the torsion

        or just get the torsion names and keys?
        build a dictionary
        (build full dct for all tors in rotor
         split dcts into subdcts for groupings, including single torsion)
    """

    tors_axes = all_axes(gra, lin_keys)

    tors_names = tuple(automol.zmat.torsion_coordinate_name(zma, *keys)
                       for keys in tors_axes)

    return dict(zip(tors_names, tors_axes))


def all_axes(gra, lin_keys):
    """ Build the torsion axes
    """
    return automol.graph.rotational_bond_keys(gra, lin_keys=lin_keys)


def all_groups(gra, lin_keys):
    """ Generate torsion groups make generalizable to multiple axes
    """

    axes = all_axes(gra, lin_keys)
    grps = ()
    for axis in axes:
        grps += (groups(gra, axis),)

    return grps


def all_symmetries(gra, lin_keys):
    """ Generate torsion groups make generalizable to multiple axes
    """

    axes = all_axes(gra, lin_keys)
    syms = ()
    for axis in axes:
        syms += (symmetry(gra, axis, lin_keys),)

    return syms


def groups(gra, axis):
    """ Generate torsion groups make generalizable to multiple axes
    """
    return automol.graph.rotational_groups(gra, *axis)


def symmetry(gra, axis, lin_keys):
    """ Obtain the symmetry number for the torsion
    """
    return automol.graph.rotational_symmetry_number(
        gra, *axis, lin_keys=lin_keys)


# I/O
def string(torsion):
    """ write the torsions to a string
    """

    def _encode_group(group_keys):
        group_str = '-'.join(str(val) for val in group_keys)
        return group_str

    tors_dct = {
        'axis': _encode_group(torsion.axis),
        'groups': _encode_group(torsion.groups),
        'symmetry': torsion.symmetry,
        # 'span': round(dct.get('span', None) * phycon.RAD2DEG, 2)
    }

    tors_str = yaml.dump(tors_dct, sort_keys=False)

    return tors_str


def from_string(tors_str):
    """ read the transformation from a string
    """

    def _decode_group(group_str):
        group_keys = tuple(map(int, group_str.split('-')))
        return group_keys

    torsions = ()

    tors_dct = yaml.load(tors_str, Loader=yaml.FullLoader)
    for name, dct in tors_dct.values():
        axis = _decode_group(dct['axis'])
        grps = _decode_group(dct['groups'])
        symm = dct

        torsions += (Torsion('', name, axis, grps, symm),)

    return torsions


def _sort_tors_names(tors_names):
    """ sort torsional names
    """
    tors_names = list(tors_names)
    tors_names.sort(key=lambda x: int(x.split('D')[1]))
    return tors_names
