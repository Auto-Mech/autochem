"""Implements a data structure for encoding rotor information

Rotors:
    zmatrix: The z-matrix relative to which the rotors are defined
    rotor_list: The list of rotors

Rotor = List[Torsion]  # each rotor is a list of torsion objects

Torsion:
    name: z-matrix coordinate name
    axis: z-matrix keys identifying the rotational axis
    groups: z-matrix keys identifying the rotational groups
    symetry: the rotational symmetry of the torsion
"""
import dataclasses
import itertools
from typing import List, Tuple, Union

import numpy
import yaml
from phydat import phycon

from automol import graph, zmat
from automol.util import zmat_conv

Axis = Tuple[int, int]
Groups = Tuple[List[int], List[int]]
Grid = List[float]


@dataclasses.dataclass
class Torsion:
    """Encodes information for a single torsion, which is one component of a rotor

    :param name: The z-matrix coordinate name
    :type name: str
    :param axis: The pair of atom keys defining the rotational axis
    :type axis: Tuple[int, int]
    :param groups: The sets of atoms keys defining the rotational groups
    :type groups: Tuple[List[int], List[int]]
    :type symmetry: The rotational symmetry number of the torsion
    :type symmetry: int
    """

    name: str
    axis: Axis
    groups: Groups
    symmetry: int


Rotor = List[Torsion]


@dataclasses.dataclass
class Rotors:
    """Encodes information for all rotors of an equilibrium or transition state
    structure

    :param zmatrix: The z-matrix
    :type zmatrix: automol zmat data structure
    :param rotors: The list of rotors for this structure
    :type rotors: List[Rotor]
    """

    zmatrix: tuple  # automol zmat data structure
    rotor_list: List[Rotor]


# Constructors
def from_data(
    zma,
    tors_lst: List[Torsion],
    tors_names_lst: List[List[str]] = None,
    multi: bool = False,
) -> Rotors:
    """Construct rotors object from existing data

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param tors_lst: A list of torsions
    :type tors_lst: List[Torsion]
    :param tors_names_lst: A list of lists of torsion names identifying which torsions
        should be included and, for multi=True, how they should be grouped;
        the dimensionality of 5+-dimensional rotors will be automatically reduced
    :type tors_names_lst: List[List[str]]
    :param multi: Construct multi-dimensional rotors? defaults to False
    :type multi: bool, optional
    :param gra: Optionally, specify the z-matrix connectivity with a graph
    :type gra: automol graph data structure
    :returns: The rotors for this structure
    :rtype: Rotors
    """

    if tors_names_lst is None:
        tors_names = [t.name for t in tors_lst]
        tors_names_lst = [tors_names] if multi else [[n] for n in tors_names]

    tors_dct = {t.name: t for t in tors_lst}

    rotor_lst = []
    for names in tors_names_lst:
        rotor = [tors_dct[n] for n in names]
        if len(rotor) > 4:
            rotor_lst.extend(_parition_high_dimensional_rotor(zma, rotor))
        else:
            rotor_lst.append(rotor)

    rotor_lst = tuple(map(tuple, rotor_lst))
    return Rotors(zmatrix=zma, rotor_list=rotor_lst)


def from_zmatrix(
    zma, gra=None, tors_names_lst: List[List[str]] = None, multi: bool = False
) -> Rotors:
    """Construct rotors by inferring torsions from the z-matrix and its graph

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param gra: Optionally, specify the z-matrix connectivity with a graph
    :type gra: automol graph data structure
    :param tors_names_lst: A list of lists of torsion names identifying which torsions
        should be included and (for multi=True) how they should be grouped;
        dimensionality of 5+-dimensional rotors will be reduced
    :type tors_names_lst: List[List[str]]
    :param multi: Construct multi-dimensional rotors? defaults to False
    :type multi: bool, optional
    :returns: The rotors for this structure
    :rtype: Rotors
    """
    gra = zmat.graph(zma, stereo=True, dummy=True) if gra is None else gra
    lin_keys = graph.linear_atom_keys(gra, dummy=True)

    rot_bkeys = graph.rotational_bond_keys(gra, lin_keys=lin_keys)

    # Read in the torsion coordinate names and sort them in z-matrix order
    tors_names = [zmat.torsion_coordinate_name(zma, *bk) for bk in rot_bkeys]
    tors_names = [n for n in zmat.dihedral_angle_names(zma) if n in tors_names]

    tors_lst = []
    for name in tors_names:
        axis = zmat.torsion_axis(zma, name)
        tors = Torsion(
            name=name,
            axis=axis,
            groups=graph.rotational_groups(gra, *axis),
            symmetry=graph.rotational_symmetry_number(gra, *axis, lin_keys=lin_keys),
        )
        tors_lst.append(tors)

    return from_data(
        zma=zma, tors_lst=tors_lst, tors_names_lst=tors_names_lst, multi=multi
    )


# Getters
def torsion_names(
    rotors: Rotors, flat: bool = False
) -> Union[List[str], List[List[str]]]:
    """Get the torsion coordinate names for a Rotors object

    :param rotors: A rotors object
    :type rotors: Rotors
    :param flat: Return a flat list instead?, defaults to False
    :type flat: bool, optional
    :return: A flat or 2D list of names
    :rtype: Union[List[str], List[List[str]]]
    """
    names_lst = [[t.name for t in r] for r in rotors.rotor_list]
    names = tuple(itertools.chain(*names_lst)) if flat else tuple(map(tuple, names_lst))
    return names


def torsion_axes(
    rotors: Rotors, key_typ: str = "zmat", flat: bool = False
) -> Union[List[Axis], List[List[Axis]]]:
    """Get the rotational axes for a Rotors object

    :param rotors: A rotors object
    :type rotors: Rotors
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param flat: Return a flat list instead?, defaults to False
    :type flat: bool, optional
    :return: A flat or 2D list of rotational axes
    :rtype: Union[List[Axis], List[List[Axis]]]
    """
    assert key_typ in ("geom", "zmat"), f"Invalid key type {key_typ} requested"

    axes_lst = [[t.axis for t in r] for r in rotors.rotor_list]

    if key_typ == "geom":
        zc_ = zmat.conversion_info(rotors.zmatrix)
        axes_lst = zmat_conv.geometry_keys(zc_, axes_lst)

    axes = tuple(itertools.chain(*axes_lst)) if flat else tuple(map(tuple, axes_lst))
    return axes


def torsion_groups(
    rotors: Rotors, key_typ: str = "zmat", flat: bool = False
) -> Union[List[Groups], List[List[Groups]]]:
    """Get the rotational groups for a Rotors object

    :param rotors: A rotors object
    :type rotors: Rotors
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param flat: Return a flat list instead?, defaults to False
    :type flat: bool, optional
    :return: A flat or 2D list of rotational groups
    :rtype: Union[List[Groups], List[List[Groups]]]
    """
    assert key_typ in ("geom", "zmat"), f"Invalid key type {key_typ} requested"

    groups_lst = [[t.groups for t in r] for r in rotors.rotor_list]

    if key_typ == "geom":
        zc_ = zmat.conversion_info(rotors.zmatrix)
        groups_lst = zmat_conv.geometry_keys(zc_, groups_lst)

    groups = (
        tuple(itertools.chain(*groups_lst)) if flat else tuple(map(tuple, groups_lst))
    )
    return groups


def torsion_symmetries(
    rotors: Rotors, flat: bool = False
) -> Union[List[int], List[List[int]]]:
    """Get the rotational symmetries for a Rotors object

    :param rotors: A rotors object
    :type rotors: Rotors
    :param flat: Return a flat list instead?, defaults to False
    :type flat: bool, optional
    :return: A flat or 2D list of rotational symmetries
    :rtype: Union[List[int], List[List[int]]]
    """
    syms = tuple(tuple(t.symmetry for t in r) for r in rotors.rotor_list)
    if flat:
        syms = tuple(itertools.chain(*syms))
    return syms


def torsion_grids(
    rotors: Rotors, span=2 * numpy.pi, increment=30 * phycon.DEG2RAD, flat: bool = False
) -> Union[List[Grid], List[List[Grid]]]:
    """Get the rotational grids for a Rotors object

    :param rotors: A rotors object
    :type rotors: Rotors
    :param flat: Return a flat list instead?, defaults to False
    :type flat: bool, optional
    :return: A flat or 2D list of rotational grids
    :rtype: Union[List[Grid], List[List[Grid]]]
    """
    zma = rotors.zmatrix

    def grid_(tors: Torsion) -> Grid:
        symm = tors.symmetry
        # [0, 30, 60, 90, ...] << in degrees
        grid = numpy.arange(0, span / symm, increment)
        # Start from the equilibrium value
        grid += zmat.value(zma, tors.name)
        return tuple(map(float, grid))

    grids = tuple(tuple(map(grid_, r)) for r in rotors.rotor_list)
    if flat:
        grids = tuple(itertools.chain(*grids))
    return grids


def torsion_dimensions(rotors: Rotors) -> List[int]:
    """Get the dimensions for the rotors in a Rotors object

    :param rotors: A rotors object
    :type rotors: Rotors
    :return: The dimension of each rotor in order
    :rtype: List[int]
    """
    return tuple(map(len, rotors.rotor_list))


def torsion_list(rotors: Rotors) -> List[Torsion]:
    """Get the torsions in a rotors object as a flat list

    :param rotors: A rotors object
    :type rotors: Rotors
    :return: The torsions in the rotors object, sorted in z-matrix order
    :rtype: List[Torsion]
    """
    tors_dct = {t.name: t for r in rotors.rotor_list for t in r}
    tors_lst = tuple(
        tors_dct[n] for n in zmat.dihedral_angle_names(rotors.zmatrix) if n in tors_dct
    )
    return tors_lst


# Torsion list functions
def torsion_list_string(tors_lst: List[Torsion], one_indexed: bool = True) -> str:
    """Write a list of torsions to a string

    :param tors_lst: A list of torsions
    :type tors_lst: List[Torsion]
    :param one_indexed: Is this a one-indexed string? defaults to True
    :type one_indexed: bool, optional
    :returns: A string representations of the torsions, as a flattened list
    :rtype: str
    """
    tors_yml_dct = torsion_list_yaml_data(tors_lst, one_indexed=one_indexed)
    tors_str = yaml.dump(tors_yml_dct, sort_keys=False)
    return tors_str


def torsion_list_from_string(tors_str: str, one_indexed: bool = True) -> List[Torsion]:
    """Write a list of torsions to a string

    :param tors_str: A string representations of the torsions, as a flattened list
    :type tors_str: str
    :param one_indexed: Is this a one-indexed string? defaults to True
    :type one_indexed: bool, optional
    :returns: A list of torsions
    :rtype: List[Torsion]
    """
    tors_yml_dct = yaml.load(tors_str, Loader=yaml.FullLoader)
    tors_lst = torsion_list_from_yaml_data(tors_yml_dct, one_indexed=one_indexed)
    return tors_lst


def torsion_list_yaml_data(tors_lst: List[Torsion], one_indexed: bool = True) -> dict:
    """Write a list of torsions to a yaml-formatted dictionary

    :param tors_lst: A list of torsions
    :type tors_lst: List[Torsion]
    :param one_indexed: Is this a one-indexed string? defaults to True
    :type one_indexed: bool, optional
    :returns: A string representations of the torsions, as a flattened list
    :rtype: str
    """
    shift = 1 if one_indexed else 0

    tors_yml_dct = {}
    for tors in tors_lst:
        axis1, axis2 = (k + shift for k in tors.axis)
        groups = ([k + shift for k in g] for g in tors.groups)
        group1, group2 = ("-".join(map(str, g)) if len(g) > 1 else g[0] for g in groups)
        tors_yml_dct[tors.name] = {
            "axis1": axis1,
            "group1": group1,
            "axis2": axis2,
            "group2": group2,
            "symmetry": tors.symmetry,
        }
    return tors_yml_dct


def torsion_list_from_yaml_data(tors_yml_dct: dict, one_indexed: bool = True) -> dict:
    """Read a list of torsions out of a yaml-formatted torsion dictionary

    :param tors_lst: A list of torsions
    :type tors_lst: List[Torsion]
    :param one_indexed: Is this a one-indexed string? defaults to True
    :type one_indexed: bool, optional
    :returns: A string representations of the torsions, as a flattened list
    :rtype: str
    """
    shift = -1 if one_indexed else 0

    tors_lst = []
    for name, vals_dct in tors_yml_dct.items():
        raw_axis = list(map(vals_dct.__getitem__, ["axis1", "axis2"]))
        raw_groups = list(map(vals_dct.__getitem__, ["group1", "group2"]))
        raw_groups = [
            [g] if isinstance(g, int) else map(int, g.split("-")) for g in raw_groups
        ]

        tors = Torsion(
            name=name,
            axis=tuple(k + shift for k in raw_axis),
            groups=tuple(tuple(k + shift for k in g) for g in raw_groups),
            symmetry=vals_dct["symmetry"],
        )

        tors_lst.append(tors)

    return tuple(tors_lst)


# helper functions
def _parition_high_dimensional_rotor(zma, rotor: Rotor) -> List[Rotor]:
    """Partition a high-dimensional rotor into rotors of smaller dimension"""
    # Separate methyl rotors from non-methyl rotors
    ch3_tors_lst = []
    rem_tors_lst = []
    for tors in rotor:
        if _is_methyl_rotor(zma, tors):
            ch3_tors_lst.append(tors)
        else:
            rem_tors_lst.append(tors)

    # If there are <= 4 non-methyl rotors, group them into a single, multi-dimensional
    # rotor, followed by individual methyl rotors
    if len(rem_tors_lst) <= 4:
        return [rem_tors_lst] + [[t] for t in ch3_tors_lst]

    # If removing methyl rotors didn't reduce the dimensionality down to <= 4, flatten
    # this down to a list of one-dimensional rotors
    return [[t] for t in rotor]


def _is_methyl_rotor(zma, tors: Torsion) -> bool:
    """From a z-matrix and a torsion object, identify this is a methyl rotor"""
    symbs = zmat.symbols(zma)
    group_keys = [[k] + list(g) for k, g in zip(tors.axis, tors.groups)]
    group_symbs = [list(map(symbs.__getitem__, ks)) for ks in group_keys]
    return bool(any(g == ["C", "H", "H", "H"] for g in group_symbs))
