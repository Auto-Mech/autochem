"""
 Build the rotor and torsion objects

 Rotor: (tors_obj_1, tors_obj_2, tors_obj_3)
"""

# import numpy
import yaml
import automol
# from automol.rotor._par import TorsionParam


class Torsion:
    """ Describes a torsion, which one or more make up a rotor
    """

    def __init__(self, zma, name, axis, groups, symmetry):
        """ constructor
        """

        self.name = name
        self.zma = zma
        self.symmetry = symmetry
        self.groups = groups
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


def torsion_lst(zma, gra, lin_keys):
    """  Build a list of torsion objects
    """

    # Build the torsion objects
    _name_axis_dct = name_axis_dct(zma, gra, lin_keys)

    tors_obj_lst = ()
    for name, axis in _name_axis_dct.items():
        grps = torsion_groups(gra, axis)
        symm = torsion_symmetry(gra, axis, lin_keys)
        tors_obj_lst += (Torsion(zma, name, axis, grps, symm),)

    return tors_obj_lst


def name_axis_dct(zma, gra, lin_keys):
    """ Generate the bond keys for the torsion

        or just get the torsion names and keys?
        build a dictionary
        (build full dct for all tors in rotor
         split dcts into subdcts for groupings, including single torsion)
    """

    tors_axes = all_torsion_axes(gra, lin_keys)

    tors_names = tuple(automol.zmat.torsion_coordinate_name(zma, *keys)
                       for keys in tors_axes)

    return dict(zip(tors_names, tors_axes))


def all_torsion_axes(gra, lin_keys):
    """ Build the torsion axes
    """
    return automol.graph.rotational_bond_keys(gra, lin_keys=lin_keys)


def all_torsion_groups(gra, lin_keys):
    """ Generate torsion groups make generalizable to multiple axes
    """

    axes = all_torsion_axes(gra, lin_keys)
    grps = ()
    for axis in axes:
        grps += (torsion_groups(gra, axis),)

    return grps


def all_torsion_symmetries(gra, lin_keys):
    """ Generate torsion groups make generalizable to multiple axes
    """

    axes = all_torsion_axes(gra, lin_keys)
    syms = ()
    for axis in axes:
        syms += (torsion_symmetry(gra, axis, lin_keys),)

    return syms


def torsion_groups(gra, axis):
    """ Generate torsion groups make generalizable to multiple axes
    """
    return automol.graph.rotational_groups(gra, *axis)


def torsion_symmetry(gra, axis, lin_keys):
    """ Obtain the symmetry number for the torsion
    """
    return automol.graph.rotational_symmetry_number(
        gra, *axis, lin_keys=lin_keys)


# I/O
def string(tors_lst):
    """ Write a list torsions to a string
    """

    def _encode_idxs(idxs):
        return '-'.join(str(val) for val in idxs)

    tors_dct = {}
    for tors in tors_lst:
        _axis = tors.axis
        _grps = tors.groups
        tors_dct[tors.name] = {
                'axis1': _encode_idxs(_axis[0]),
                'group1': _encode_idxs(_grps[0]),
                'axis2': _encode_idxs(_axis[1]),
                'group2': _encode_idxs(_grps[1]),
                'symmetry': tors.symmetry,
        }

    tors_str = yaml.dump(tors_dct, sort_keys=False)

    return tors_str


def from_string(tors_str):
    """ read the transformation from a string
    """

    def _decode_idxs(idxs_str):
        return tuple(map(int, idxs_str.split('-')))

    inf_dct = {}

    tors_dct = yaml.load(tors_str, Loader=yaml.FullLoader)
    for name, dct in tors_dct.items():
        _axis = (_decode_idxs(dct['axis1']), _decode_idxs(dct['axis2']))
        _grps = (_decode_idxs(dct['group1']), _decode_idxs(dct['group2']))
        symm = dct['symmetry']

        inf_dct[name] = {'axis': _axis, 'groups': _grps, 'symmetry': symm}
        # torsions += (Torsion('', name, _axis, _grps, symm),)

    return inf_dct


def _sort_tors_names(tors_names):
    """ sort torsional names
    """
    tors_names = list(tors_names)
    tors_names.sort(key=lambda x: int(x.split('D')[1]))
    return tors_names
