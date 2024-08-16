"""Implements a data structure for encoding rotor information

Rotors contain of one or more Torsion data structures
"""

import dataclasses
import itertools
from typing import Dict, List, Optional, Union

from phydat import phycon

from .. import graph, zmat
from ..util import ZmatConv
from . import potent, tors
from .potent import Potential
from .tors import Axis, DihCoord, Grid, Groups, Torsion


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


def from_data(zma, tor_lst: List[Torsion], pot: Optional[Potential] = None) -> Rotor:
    """Construct a rotor from data

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param tor_lst: A list of torsions
    :type tor_lst: List[Torsion]
    :param pot: Optionally, specify a potential for the rotor, defaults to None
    :type pot: Potential, optional
    :return: A rotor data structure
    :rtype: Rotor
    """
    zma_names = zmat.dihedral_angle_names(zma)
    tor_names = tuple(map(tors.name, tor_lst))

    assert all(
        n in zma_names for n in tor_names
    ), f"Torsion names don't match z-matrix:\n{tor_names}\n{zma}"

    if pot is not None:
        pot_names = potent.coordinate_names(pot)
        assert (
            tor_names == pot_names
        ), f"Torsion names don't match potential:\n{tor_names}\n{pot_names}"

    return Rotor(zmatrix=zma, torsions=tor_lst, potential=pot)


# Getters
def zmatrix(rotor: Rotor):
    """Get the z-matrix associated with a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :return: The z-matrix
    :rtype: automol zmat data structure
    """
    return rotor.zmatrix


def torsions(
    rotor: Rotor, key_typ: str = "zmat", zc_: Optional[ZmatConv] = None
) -> List[Torsion]:
    """Get the torsions making up a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :type zc_: Optional[ZmatConv], optional
    :return: The list of torsions
    :rtype: List[Torsion]
    """
    assert key_typ in ("geom", "zmat"), f"Invalid key type {key_typ} requested"

    tor_lst = rotor.torsions

    if key_typ == "geom":
        zc_ = zmat.conversion_info(zmatrix(rotor)) if zc_ is None else zc_
        tor_lst = tuple(tors.with_geometry_indices(t, zc_=zc_) for t in tor_lst)

    return tor_lst


def torsion_dict(
    rotor: Rotor, key_typ: str = "zmat", zc_: Optional[ZmatConv] = None
) -> Dict[str, Torsion]:
    """Get the torsions in a rotor as a dictionary, by coordinate name

    :param rotor: A rotor
    :type rotor: Rotor
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :type zc_: Optional[ZmatConv], optional
    :return: The torsions in the rotor object, by coordinate name
    :rtype: Dict[str, Torsion]
    """
    return {tors.name(t): t for t in torsions(rotor, key_typ=key_typ, zc_=zc_)}


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
    return tuple(map(tors.name, torsions(rotor)))


def torsion_coordinate(
    rotor: Rotor, key_typ: str = "zmat", zc_: Optional[ZmatConv] = None
) -> List[DihCoord]:
    """Get the list of dihedral coordinate keys for the torsions in a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :type zc_: Optional[ZmatConv], optional
    :return: The dihedral coordinates, in order
    :rtype: List[DihKey]
    """
    return tuple(map(tors.coordinate, torsions(rotor, key_typ=key_typ, zc_=zc_)))


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
    return tuple(map(tors.axis, torsions(rotor, key_typ=key_typ, zc_=zc_)))


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
    return tuple(map(tors.groups, torsions(rotor, key_typ=key_typ, zc_=zc_)))


def torsion_symmetries(rotor: Rotor) -> List[int]:
    """Get the list of rotational symmetries for the torsions in a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :return: The torsion rotational symmetries, in order
    :rtype: List[int]
    """
    return tuple(map(tors.symmetry, torsions(rotor)))


def torsion_grids(rotor: Rotor, increment: float = 30 * phycon.DEG2RAD) -> List[Grid]:
    """Get the coordinate grids for the torsions in a rotor

    :param rotor: A rotor
    :type rotor: Rotor
    :param increment: The grid increment, in radians
    :type increment: float
    :return: The torsion coordinate grids, in order
    :rtype: List[Grid]
    """
    zma = zmatrix(rotor)
    return tuple(tors.grid(t, zma=zma, increment=increment) for t in torsions(rotor))


# Setters
def set_potential(rotor: Rotor, pot: Potential, in_place: bool = False) -> Rotor:
    """Set the rotor potential

    :param rotor: A rotor
    :type rotor: Rotor
    :param pot: The rotor potential
    :type pot: Potential
    :param in_place: Do this in-place, mutating the object? defaults to False
    :type in_place: bool, optional
    :return: A new rotor
    :rtype: Rotor
    """
    if in_place:
        rotor.potential = pot
    else:
        rotor = from_data(zma=zmatrix(rotor), tor_lst=torsions(rotor), pot=pot)
    return rotor


# Transformations
def partition_high_dimensional_rotor(rotor: Rotor) -> List[Rotor]:
    """Partition a high-dimensional rotor into rotors of smaller dimension"""
    zma = zmatrix(rotor)
    symbs = zmat.symbols(zma)

    def is_methyl_rotor_(tor: Torsion) -> bool:
        """From a z-matrix and a torsion object, identify if this is a methyl rotor"""
        group_keys = [[k] + list(g) for k, g in zip(tors.axis(tor), tors.groups(tor))]
        group_symbs = [list(map(symbs.__getitem__, ks)) for ks in group_keys]
        return bool(any(g == ["C", "H", "H", "H"] for g in group_symbs))

    tor_lst = torsions(rotor)

    if len(tor_lst) <= 4:
        return (rotor,)

    # Separate methyl rotors from non-methyl rotors
    ch3_lst = []
    rem_lst = []
    for tor in tor_lst:
        if is_methyl_rotor_(tor):
            ch3_lst.append(tor)
        else:
            rem_lst.append(tor)

    if len(rem_lst) <= 4:
        # If there are <= 4 non-methyl rotors, group them into a single,
        # multi-dimensional rotor, followed by individual methyl rotors
        tor_lsts = [rem_lst] + [[t] for t in ch3_lst]
    else:
        # If removing methyl rotors didn't reduce the dimensionality down to <= 4,
        # flatten this down to a list of one-dimensional rotors
        tor_lsts = [[t] for t in tor_lst]

    rotors = tuple(from_data(zma=zma, tor_lst=ts) for ts in tor_lsts)
    return rotors


# Rotor List functions
# # Constructors
def rotors_from_data(
    zma,
    tor_lst: List[Torsion],
    tor_names_lst: List[List[str]] = None,
    partition: bool = True,
    multi: bool = False,
) -> List[Rotor]:
    """Construct rotors object from existing data

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param tor_lst: A list of torsions
    :type tor_lst: List[Torsion]
    :param tor_names_lst: A list of lists of torsion names identifying which torsions
        should be included and how they should be grouped
    :type tor_names_lst: List[List[str]]
    :param partition: Partition high-dimensional rotors? defaults to True
    :type partition: bool, optional
    :param multi: If no grouping was specified, group all torsions into a single rotor?
        Will be partitioned to reduce dimensionality if `partition=True` is set
    :type multi: bool, optional
    :returns: A list of rotors for this z-matrix
    :rtype: List[Rotor]
    """
    if tor_names_lst is None:
        tor_names = list(map(tors.name, tor_lst))
        tor_names_lst = [tor_names] if multi else [[n] for n in tor_names]

    tor_dct = {tors.name(t): tors.update_against_zmatrix(t, zma) for t in tor_lst}

    rotors = []
    for names in tor_names_lst:
        rotor = Rotor(zmatrix=zma, torsions=[tor_dct[n] for n in names])
        rotors.append(rotor)

    if partition:
        rotors_iter = map(partition_high_dimensional_rotor, rotors)
        rotors = tuple(itertools.chain(*rotors_iter))

    return tuple(rotors)


def rotors_from_zmatrix(
    zma,
    gra=None,
    tor_names_lst: List[List[str]] = None,
    partition: bool = True,
    multi: bool = False,
) -> List[Rotor]:
    """Construct rotors by inferring torsions from the z-matrix and its graph

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param gra: Optionally, specify the z-matrix connectivity with a graph
    :type gra: automol graph data structure
    :param tor_names_lst: A list of lists of torsion names identifying which torsions
        should be included and how they should be grouped
    :type tor_names_lst: List[List[str]]
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
    znkeys_dct = zmat.neighbor_keys(zma)

    rot_bkeys = graph.rotational_bond_keys(gra, lin_keys=lin_keys)

    # Read in the torsion coordinate names and sort them in z-matrix order
    tor_names = [zmat.torsion_coordinate_name(zma, *bk) for bk in rot_bkeys]
    tor_names = [n for n in zmat.dihedral_angle_names(zma) if n in tor_names]

    tor_lst = []
    for name in tor_names:
        tor_keys = list(reversed(zmat.coordinate(zma, name)))
        tor_axis = tor_keys[1:3]
        tor = tors.from_data(
            name_=name,
            coo=tor_keys,
            grps=graph.rotational_groups(gra, *tor_axis),
            symm=graph.rotational_symmetry_number(gra, *tor_axis, lin_keys=lin_keys),
            ngrps=list(map(znkeys_dct.get, tor_axis)),
        )
        tor_lst.append(tor)

    return rotors_from_data(
        zma=zma,
        tor_lst=tor_lst,
        tor_names_lst=tor_names_lst,
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
    if not rotors:
        return None

    zmas = list(map(zmatrix, rotors))
    zma, *zmas = zmas
    assert all(zma == z for z in zmas)
    return zma


def rotors_torsion_dict(
    rotors: List[Rotor], key_typ: str = "zmat"
) -> Dict[str, Torsion]:
    """Get the torsions in a list of rotors as a dictionary, by coordinate name

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :return: The torsions in the rotor object, by coordinate name
    :rtype: Dict[str, Torsion]
    """
    if not rotors:
        return {}

    zc_ = None if key_typ == "zmat" else zmat.conversion_info(rotors_zmatrix(rotors))
    return {
        tors.name(t): t for r in rotors for t in torsions(r, key_typ=key_typ, zc_=zc_)
    }


def rotors_torsions(
    rotors: List[Rotor], key_typ: str = "zmat", flat: bool = False, sort: bool = False
) -> Union[List[Torsion], List[List[Torsion]]]:
    """Get the torsion coordinate names from a list of rotors

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param flat: Return a flat list instead of grouping by rotors?, defaults to False
    :type flat: bool, optional
    :param sort: Return a flat list sorted in z-matrix order?, defaults to False
    :type sort: bool, optional
    :return: A flat or grouped list of torsion names
    :rtype: Union[List[str], List[List[str]]]
    """
    if not rotors:
        return ()

    if sort:
        tor_dct = rotors_torsion_dict(rotors, key_typ=key_typ)
        zma = rotors_zmatrix(rotors)
        return tuple(tor_dct[n] for n in zmat.dihedral_angle_names(zma) if n in tor_dct)

    zc_ = zmat.conversion_info(rotors_zmatrix(rotors)) if key_typ == "geom" else None
    tor_lst = tuple(torsions(r, key_typ=key_typ, zc_=zc_) for r in rotors)
    return tuple(itertools.chain(*tor_lst)) if flat else tor_lst


def rotors_potentials(rotors: List[Rotor]) -> List[Optional[Potential]]:
    """Get the rotor potentials from a list of rotors, in order

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :return: The potentials, in order
    :rtype: List[Optional[Potential]]
    """
    return tuple(map(potential, rotors))


def rotors_have_potentials(rotors: List[Rotor]) -> bool:
    """Do these rotors have potentials?

    :param rotors: A list of rotor objects
    :return: `True` if they do, `False` if they don't
    """
    return all(map(potent.has_defined_values, rotors_potentials(rotors)))


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


def rotors_torsion_coordinates(
    rotors: List[Rotor], key_typ: str = "zmat", flat: bool = False
) -> Union[List[DihCoord], List[List[DihCoord]]]:
    """Get the torsion rotational keys from a list of rotors

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :param flat: Return a flat list instead of grouping by rotors?, defaults to False
    :type flat: bool, optional
    :return: A flat or grouped list of torsion keys
    :rtype: Union[List[DihKey], List[List[DihKey]]]
    """
    zc_ = None if key_typ == "zmat" else zmat.conversion_info(rotors_zmatrix(rotors))
    coo_lst = tuple(torsion_coordinate(r, key_typ=key_typ, zc_=zc_) for r in rotors)
    return tuple(itertools.chain(*coo_lst)) if flat else coo_lst


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
    rotors: List[Rotor], flat: bool = False, increment: float = 30 * phycon.DEG2RAD
) -> Union[List[Grid], List[List[Grid]]]:
    """Get the torsion coordinate grids from a list of rotors

    :param rotors: A list of rotor objects
    :type rotors: List[Rotor]
    :param flat: Return a flat list instead of grouping by rotors?, defaults to False
    :type flat: bool, optional
    :param increment: The grid increment, in radians
    :type increment: float
    :return: A flat or grouped list of torsion grids
    :rtype: Union[List[Grid], List[List[Grid]]]
    """
    grids_lst = [torsion_grids(r, increment=increment) for r in rotors]
    return tuple(itertools.chain(*grids_lst)) if flat else grids_lst
