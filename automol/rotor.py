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
from typing import Dict, List, Optional, Tuple, Union

import numpy
import yaml
from phydat import phycon

from automol import graph, potent, zmat
from automol.potent import Potential
from automol.util import ZmatConv, zmat_conv

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


@dataclasses.dataclass
class Rotor:
    """Encodes information for a single rotor, consisting of a set of torsions
    associated with a potential and a z-matrix

    :param zmatrix: The z-matrix
    :type zmatrix: automol zmat data structure
    :param torsions: The torsions constituting this rotor
    :type torsions: List[Torsion]
    :param potential: The rotor potential, defaults to None
    :type potential: Optional[Potential]
    """

    zmatrix: tuple  # automol zmat data structure
    torsions: List[Torsion]
    potential: Optional[Potential] = None


def from_data(zma, tors_lst: List[Torsion], pot: Optional[Potential] = None) -> Rotor:
    """Construct a rotor from data

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param tors_lst: A list of torsions
    :type tors_lst: List[Torsion]
    :param pot: Optionally, specify a potential for the rotor, defaults to None
    :type pot: Potential, optional
    :return: A rotor data structure
    :rtype: Rotor
    """
    zma_names = zmat.dihedral_angle_names(zma)
    tors_names = tuple(t.name for t in tors_lst)

    assert all(
        n in zma_names for n in tors_names
    ), f"Torsion names don't match z-matrix:\n{tors_names}\n{zma}"

    if pot is not None:
        pot_names = potent.coordinate_names(pot)
        assert (
            tors_names == pot_names
        ), f"Torsion names don't match potential:\n{tors_names}\n{pot_names}"

    return Rotor(zmatrix=zma, torsions=tors_lst, potential=pot)


# Getters
def zmatrix(rotor: Rotor):
    """Get the z-matrix associated with a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :return: The z-matrix
    :rtype: automol zmat data structure
    """
    return rotor.zmatrix


def torsions(rotor: Rotor) -> List[Torsion]:
    """Get the torsions making up a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :return: The list of torsions
    :rtype: List[Torsion]
    """
    return rotor.torsions


def torsion_dict(rotor: Rotor) -> Dict[str, Torsion]:
    """Get the torsions in a rotor as a dictionary, by coordinate name

    :param rotor: A rotor
    :type rotor: Rotor
    :return: The torsions in the rotor object, by coordinate name
    :rtype: Dict[str, Torsion]
    """
    return {t.name: t for t in torsions(rotor)}


def potential(rotor: Rotor) -> Optional[Potential]:
    """Get the rotor potential, if set

    :param rotor: A rotor
    :type rotor: Rotor
    :return: The rotor potential, if set, otherwise None
    :rtype: Optional[Potential]
    """
    return rotor.potential


def dimension(rotor: Rotor) -> int:
    """Get the dimension of (number of torsions in) a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :return: The dimension
    :rtype: int
    """
    return len(torsions(rotor))


def torsion_names(rotor: Rotor) -> List[str]:
    """Get the list of coordinate names for the torsions in a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :return: The torsion coordinate names, in order
    :rtype: List[str]
    """
    return tuple(t.name for t in torsions(rotor))


def torsion_axes(
    rotor: Rotor, key_typ: str = "zmat", zc_: Optional[ZmatConv] = None
) -> List[Axis]:
    """Get the list of rotational axes for the torsions in a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :type zc_: Optional[ZmatConv], optional
    :return: The torsion rotational axes, in order
    :rtype: List[Axis]
    """
    assert key_typ in ("geom", "zmat"), f"Invalid key type {key_typ} requested"

    axes = tuple(t.axis for t in torsions(rotor))

    if key_typ == "geom":
        zc_ = zmat.conversion_info(zmatrix(rotor)) if zc_ is None else zc_
        axes = zmat_conv.geometry_keys(zc_, axes)

    return axes


def torsion_groups(
    rotor: Rotor, key_typ: str = "zmat", zc_: Optional[ZmatConv] = None
) -> List[Groups]:
    """Get the list of rotational groups for the torsions in a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :type zc_: Optional[ZmatConv], optional
    :return: The torsion rotational groups, in order
    :rtype: List[Groups]
    """
    assert key_typ in ("geom", "zmat"), f"Invalid key type {key_typ} requested"

    groups_lst = tuple(t.groups for t in torsions(rotor))

    if key_typ == "geom":
        zc_ = zmat.conversion_info(zmatrix(rotor)) if zc_ is None else zc_
        groups_lst = zmat_conv.geometry_keys(zc_, groups_lst)

    return groups_lst


def torsion_symmetries(rotor: Rotor) -> List[int]:
    """Get the list of rotational symmetries for the torsions in a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :return: The torsion rotational symmetries, in order
    :rtype: List[int]
    """
    return tuple(t.symmetry for t in torsions(rotor))


def torsion_grids(
    rotor: Rotor, span=2 * numpy.pi, increment=30 * phycon.DEG2RAD
) -> List[Grid]:
    """Get the coordinate grids for the torsions in a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :return: The torsion coordinate grids, in order
    :rtype: List[Grid]
    """
    zma = zmatrix(rotor)

    def grid_(tors: Torsion) -> Grid:
        symm = tors.symmetry
        # [0, 30, 60, 90, ...] << in degrees
        grid = numpy.arange(0, span / symm, increment)
        # Start from the equilibrium value
        grid += zmat.value(zma, tors.name)
        return tuple(map(float, grid))

    return tuple(map(grid_, torsions(rotor)))


# Setters
def set_potential(rotor: Rotor, pot: Potential) -> Rotor:
    """Set the rotor potential

    :param rotor: A rotor
    :type rotor: Rotor
    :param pot: The rotor potential
    :type pot: Potential
    :return: A new rotor
    :rtype: Rotor
    """
    return from_data(zma=zmatrix(rotor), tors_lst=torsions(rotor), pot=pot)


# Transformations
def partition_high_dimensional_rotor(rotor: Rotor) -> List[Rotor]:
    """Partition a high-dimensional rotor into rotors of smaller dimension"""
    zma = zmatrix(rotor)
    symbs = zmat.symbols(zma)

    def is_methyl_rotor_(tors: Torsion) -> bool:
        """From a z-matrix and a torsion object, identify if this is a methyl rotor"""
        group_keys = [[k] + list(g) for k, g in zip(tors.axis, tors.groups)]
        group_symbs = [list(map(symbs.__getitem__, ks)) for ks in group_keys]
        return bool(any(g == ["C", "H", "H", "H"] for g in group_symbs))

    tors_lst = torsions(rotor)

    if len(tors_lst) <= 4:
        return (rotor,)

    # Separate methyl rotors from non-methyl rotors
    ch3_lst = []
    rem_lst = []
    for tors in tors_lst:
        if is_methyl_rotor_(tors):
            ch3_lst.append(tors)
        else:
            rem_lst.append(tors)

    if len(rem_lst) <= 4:
        # If there are <= 4 non-methyl rotors, group them into a single,
        # multi-dimensional rotor, followed by individual methyl rotors
        tors_lsts = [rem_lst] + [[t] for t in ch3_lst]
    else:
        # If removing methyl rotors didn't reduce the dimensionality down to <= 4,
        # flatten this down to a list of one-dimensional rotors
        tors_lsts = [[t] for t in tors_lst]

    rotors = tuple(from_data(zma=zma, tors_lst=ts) for ts in tors_lsts)
    return rotors


# Rotor List functions
# # Constructors
def rotors_from_data(
    zma,
    tors_lst: List[Torsion],
    tors_names_lst: List[List[str]] = None,
    partition: bool = True,
    multi: bool = False,
) -> List[Rotor]:
    """Construct rotors object from existing data

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param tors_lst: A list of torsions
    :type tors_lst: List[Torsion]
    :param tors_names_lst: A list of lists of torsion names identifying which torsions
        should be included and how they should be grouped
    :type tors_names_lst: List[List[str]]
    :param partition: Partition high-dimensional rotors? defaults to True
    :type partition: bool, optional
    :param multi: If no grouping was specified, group all torsions into a single rotor?
        Will be partitioned to reduce dimensionality if `partition=True` is set
    :type multi: bool, optional
    :returns: A list of rotors for this z-matrix
    :rtype: List[Rotor]
    """
    if tors_names_lst is None:
        tors_names = [t.name for t in tors_lst]
        tors_names_lst = [tors_names] if multi else [[n] for n in tors_names]

    tors_dct = {t.name: t for t in tors_lst}

    rotors = []
    for names in tors_names_lst:
        rotor = Rotor(zmatrix=zma, torsions=[tors_dct[n] for n in names])
        rotors.append(rotor)

    if partition:
        rotors_iter = map(partition_high_dimensional_rotor, rotors)
        rotors = tuple(itertools.chain(*rotors_iter))

    return tuple(rotors)


def rotors_from_zmatrix(
    zma,
    gra=None,
    tors_names_lst: List[List[str]] = None,
    partition: bool = True,
    multi: bool = False,
) -> List[Rotor]:
    """Construct rotors by inferring torsions from the z-matrix and its graph

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param gra: Optionally, specify the z-matrix connectivity with a graph
    :type gra: automol graph data structure
    :param tors_names_lst: A list of lists of torsion names identifying which torsions
        should be included and how they should be grouped
    :type tors_names_lst: List[List[str]]
    :param partition: Partition high-dimensional rotors? defaults to True
    :type partition: bool, optional
    :param multi: If no grouping was specified, group all torsions into a single rotor?
        Will be partitioned to reduce dimensionality if `partition=True` is set
    :type multi: bool, optional
    :returns: A list of rotors for this z-matrix
    :rtype: List[Rotor]
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

    return rotors_from_data(
        zma=zma,
        tors_lst=tors_lst,
        tors_names_lst=tors_names_lst,
        partition=partition,
        multi=multi,
    )


# # Getters
def rotors_zmatrix(rotors: List[Rotor]):
    """Get the z-matrix associated with a list of rotors

    Requires that all z-matrices match exactly

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :returns: The z-matrix
    :rtype: automol zmat data structure
    """
    zmas = list(map(zmatrix, rotors))
    zma, *zmas = zmas
    assert all(zma == z for z in zmas)
    return zma


def rotors_torsion_dict(rotors: List[Rotor]) -> Dict[str, Torsion]:
    """Get the torsions in a list of rotors as a dictionary, by coordinate name

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :return: The torsions in the rotor object, by coordinate name
    :rtype: Dict[str, Torsion]
    """
    return {t.name: t for r in rotors for t in torsions(r)}


def rotors_torsions(
    rotors: List[Rotor], flat: bool = False, sort: bool = False
) -> Union[List[Torsion], List[List[Torsion]]]:
    """Get the torsion coordinate names from a list of rotors

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :param flat: Return a flat list instead of grouping by rotors?, defaults to False
    :type flat: bool, optional
    :param sort: Return a flat list sorted in z-matrix order?, defaults to False
    :type sort: bool, optional
    :return: A flat or grouped list of torsion names
    :rtype: Union[List[str], List[List[str]]]
    """
    if sort:
        tors_dct = rotors_torsion_dict(rotors)
        zma = rotors_zmatrix(rotors)
        return tuple(
            tors_dct[n] for n in zmat.dihedral_angle_names(zma) if n in tors_dct
        )

    tors_lst = tuple(map(torsions, rotors))
    return tuple(itertools.chain(*tors_lst)) if flat else tors_lst


def rotors_potentials(rotors: List[Rotor]) -> List[Optional[Potential]]:
    """Get the rotor potentials from a list of rotors, in order

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :return: The potentials, in order
    :rtype: List[Optional[Potential]]
    """
    return tuple(map(potential, rotors))


def rotors_dimensions(rotors: List[Rotor]) -> List[int]:
    """Get the rotor dimensions from a list of rotors, in order

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :return: The dimensions, in order
    :rtype: List[int]
    """
    return tuple(map(dimension, rotors))


def rotors_torsion_names(
    rotors: List[Rotor], flat: bool = False
) -> Union[List[str], List[List[str]]]:
    """Get the torsion coordinate names from a list of rotors

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :param flat: Return a flat list instead of grouping by rotors?, defaults to False
    :type flat: bool, optional
    :return: A flat or grouped list of torsion names
    :rtype: Union[List[str], List[List[str]]]
    """
    names_lst = tuple(map(torsion_names, rotors))
    return tuple(itertools.chain(*names_lst)) if flat else names_lst


def rotors_torsion_axes(
    rotors: List[Rotor], key_typ: str = "zmat", flat: bool = False
) -> Union[List[Axis], List[List[Axis]]]:
    """Get the torsion rotational axes from a list of rotors

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :param flat: Return a flat list instead of grouping by rotors?, defaults to False
    :type flat: bool, optional
    :return: A flat or grouped list of torsion axes
    :rtype: Union[List[Axis], List[List[Axis]]]
    """
    zc_ = None if key_typ == "zmat" else zmat.conversion_info(rotors_zmatrix(rotors))
    axes_lst = tuple(torsion_axes(r, key_typ=key_typ, zc_=zc_) for r in rotors)
    return tuple(itertools.chain(*axes_lst)) if flat else axes_lst


def rotors_torsion_groups(
    rotors: List[Rotor], key_typ: str = "zmat", flat: bool = False
) -> Union[List[Groups], List[List[Groups]]]:
    """Get the torsion rotational groups from a list of rotors

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :param flat: Return a flat list instead of grouping by rotors?, defaults to False
    :type flat: bool, optional
    :return: A flat or grouped list of torsion groups
    :rtype: Union[List[Groups], List[List[Groups]]]
    """
    zc_ = None if key_typ == "zmat" else zmat.conversion_info(rotors_zmatrix(rotors))
    groups_lst = tuple(torsion_groups(r, key_typ=key_typ, zc_=zc_) for r in rotors)
    return tuple(itertools.chain(*groups_lst)) if flat else groups_lst


def rotors_torsion_symmetries(
    rotors: List[Rotor], flat: bool = False
) -> Union[List[str], List[List[str]]]:
    """Get the torsion coordinate symmetries from a list of rotors

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :param flat: Return a flat list instead of grouping by rotors?, defaults to False
    :type flat: bool, optional
    :return: A flat or grouped list of torsion symmetries
    :rtype: Union[List[str], List[List[str]]]
    """
    symms_lst = tuple(map(torsion_symmetries, rotors))
    return tuple(itertools.chain(*symms_lst)) if flat else symms_lst


def rotors_torsion_grids(
    rotors: List[Rotor], flat: bool = False
) -> Union[List[Grid], List[List[Grid]]]:
    """Get the torsion coordinate grids from a list of rotors

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :param flat: Return a flat list instead of grouping by rotors?, defaults to False
    :type flat: bool, optional
    :return: A flat or grouped list of torsion grids
    :rtype: Union[List[Grid], List[List[Grid]]]
    """
    grids_lst = tuple(map(torsion_grids, rotors))
    return tuple(itertools.chain(*grids_lst)) if flat else grids_lst


# Torsion List functions
def torsions_string(tors_lst: List[Torsion], one_indexed: bool = True) -> str:
    """Write a list of torsions to a string

    :param tors_lst: A list of torsions
    :type tors_lst: List[Torsion]
    :param one_indexed: Is this a one-indexed string? defaults to True
    :type one_indexed: bool, optional
    :returns: A string representations of the torsions, as a flattened list
    :rtype: str
    """
    tors_yml_dct = torsions_yaml_data(tors_lst, one_indexed=one_indexed)
    tors_str = yaml.dump(tors_yml_dct, sort_keys=False)
    return tors_str


def torsions_from_string(tors_str: str, one_indexed: bool = True) -> List[Torsion]:
    """Write a list of torsions to a string

    :param tors_str: A string representations of the torsions, as a flattened list
    :type tors_str: str
    :param one_indexed: Is this a one-indexed string? defaults to True
    :type one_indexed: bool, optional
    :returns: A list of torsions
    :rtype: List[Torsion]
    """
    tors_yml_dct = yaml.load(tors_str, Loader=yaml.FullLoader)
    tors_lst = torsions_from_yaml_data(tors_yml_dct, one_indexed=one_indexed)
    return tors_lst


def torsions_yaml_data(tors_lst: List[Torsion], one_indexed: bool = True) -> dict:
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


def torsions_from_yaml_data(tors_yml_dct: dict, one_indexed: bool = True) -> dict:
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
