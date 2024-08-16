"""Implements a data structure for encoding information about a single torsion

Torsions are used to construct Rotor data structures
"""

import dataclasses
from typing import Any

import numpy
import yaml

from phydat import phycon

from .. import zmat
from ..util import ZmatConv, zmat_conv

Axis = tuple[int, int]
DihCoord = tuple[int, int, int, int]
Group = list[int]
Groups = tuple[Group, Group]
Grid = list[float]


@dataclasses.dataclass
class Torsion:
    """Encodes information for a single torsion, which is one component of a rotor.

    :param name: The z-matrix coordinate name
    :param coordinate: The z-matrix keys defining the torsion coordinate
    :param groups: The sets of atoms keys defining the rotational groups
    :param symmetry: The rotational symmetry number of the torsion
    :param neighbor_groups: The subsets of rotational group atoms that neighbor the axis
    """

    name: str
    coordinate: DihCoord
    groups: Groups
    symmetry: int
    neighbor_groups: Groups | None = None


# Torsion functions
# # Constructors
def from_data(
    name_: str, coo: DihCoord, grps: Groups, symm: int, ngrps: Groups | None = None
) -> Torsion:
    """Construct a torsion from data.

    :param name_: The z-matrix coordinate name
    :param coo: The z-matrix keys defining the torsion coordinate
    :param groups: The sets of atoms keys defining the rotational groups
    :type symm: The rotational symmetry number of the torsion
    :return: The torsion data structure
    """
    assert len(coo) == 4, f"Invalid torsion coordinate: {coo}"
    assert len(grps) == 2, f"Invalid torsion groups: {grps}"
    return Torsion(
        name=str(name_),
        coordinate=tuple(coo),
        groups=tuple(map(tuple, grps)),
        symmetry=int(symm),
        neighbor_groups=None if ngrps is None else tuple(map(tuple, ngrps)),
    )


# # Getters
def name(tor: Torsion) -> str:
    """Get the coordinate name of a torsion.

    :param tor: A torsion
    :return: The coordinate name
    """
    return tor.name


def coordinate(
    tor: Torsion,
    key_typ: str = "zmat",
    zc_: ZmatConv | None = None,
    replace_dummy: bool = True,
) -> DihCoord:
    """Get the torsion coordinate keys.

    :param tor: A torsion
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :param replace_dummy: Replace dummy atom with other (presumably in-line) neighbor?
    :return: The torsion rotational keys
    """
    coo = tor.coordinate
    if key_typ == "geom":
        if replace_dummy:
            dkeys = zmat_conv.dummy_keys(zc_)
            ngrps = neighbor_groups(tor, dummy=False, zc_=zc_)
            end_keys0 = [coo[0], coo[-1]]
            end_keys = [
                k if k not in dkeys else next(iter(nks))
                for k, nks in zip(end_keys0, ngrps, strict=True)
            ]
            coo = (end_keys[0], coo[1], coo[2], end_keys[1])
        coo = zmat_conv.relabel_zmatrix_key_sequence(zc_, coo, dummy=True)
    return coo


def groups(tor: Torsion, key_typ: str = "zmat", zc_: ZmatConv | None = None) -> Groups:
    """Get the rotational groups of a torsion.

    :param tor: A torsion
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :return: The torsion rotational groups
    """
    grps = tor.groups
    if key_typ == "geom":
        grps = zmat_conv.relabel_zmatrix_key_sequence(zc_, grps)
    return grps


def symmetry(tor: Torsion) -> int:
    """Get the rotational symmetry of a torsion.

    :param tor: A torsion
    :return: The rotational symmetry
    """
    return tor.symmetry


def neighbor_groups(
    tor: Torsion, key_typ: str = "zmat", zc_: ZmatConv | None = None, dummy: bool = True
) -> Groups:
    """Get the neighbor groups of a torsion.

    :param tor: A torsion
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :param dummy: Include dummy atom neighbors?
    :return: The neighbor groups
    """
    ngrps = tor.neighbor_groups
    if not dummy:
        dkeys = zmat_conv.dummy_keys(zc_)
        ngrps = tuple([k for k in g if k not in dkeys] for g in ngrps)

    if key_typ == "geom":
        ngrps = zmat_conv.relabel_zmatrix_key_sequence(zc_, ngrps)
    return ngrps


# setters
def set_name(tor: Torsion, name_: str) -> Torsion:
    """Set the coordinate name of a torsion.

    :param tor: A torsion
    :param name_: The coordinate name
    :return: The torsion
    """
    return from_data(
        name_=name_,
        coo=coordinate(tor),
        grps=groups(tor),
        symm=symmetry(tor),
        ngrps=neighbor_groups(tor),
    )


def set_coordinate(tor: Torsion, coo: DihCoord) -> Torsion:
    """Set the torsion coordinate keys.

    :param tor: A torsion
    :param coo: The torsion coordinate keys
    :return: The torsion
    """
    return from_data(
        name_=name(tor),
        coo=coo,
        grps=groups(tor),
        symm=symmetry(tor),
        ngrps=neighbor_groups(tor),
    )


def set_groups(tor: Torsion, grps: Groups) -> Torsion:
    """Set the rotational groups of a torsion.

    :param tor: A torsion
    :param grps: The torsion rotational groups
    :return: The torsion
    """
    return from_data(
        name_=name(tor),
        coo=coordinate(tor),
        grps=grps,
        symm=symmetry(tor),
        ngrps=neighbor_groups(tor),
    )


def set_symmetry(tor: Torsion, symm: int) -> Torsion:
    """Set the rotational symmetry of a torsion.

    :param tor: A torsion
    :param symm: The rotational symmetry
    :return: The torsion
    """
    return from_data(
        name_=name(tor),
        coo=coordinate(tor),
        grps=groups(tor),
        symm=symm,
        ngrps=neighbor_groups(tor),
    )


def set_neighbor_groups(tor: Torsion, ngrps: Groups) -> Torsion:
    """Set the neighbor groups of a torsion.

    :param tor: A torsion
    :param ngrps: The neighbor groups
    :return: The torsion
    """
    return from_data(
        name_=name(tor),
        coo=coordinate(tor),
        grps=groups(tor),
        symm=symmetry(tor),
        ngrps=ngrps,
    )


# properties
def axis(tor: Torsion, key_typ: str = "zmat", zc_: ZmatConv | None = None) -> Axis:
    """Get the rotational axis of a torsion.

    :param tor: A torsion
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :return: The torsion rotational axis
    """
    ks_ = coordinate(tor, key_typ=key_typ, zc_=zc_)
    return ks_[1:3]


def span(tor: Torsion) -> float:
    """Get the angular span of a torsion, based on the symmetry number.

    :param tor: A torsion
    :return: The angular span of the torsion, 2 pi / symmetry
    """
    return 2 * numpy.pi / symmetry(tor)


def grid(tor: Torsion, zma, increment: float = 30 * phycon.DEG2RAD) -> Grid:
    """Get the coordinate grid for a torsion.

    :param tor: A torsion
    :param zma: The z-matrix associated with this torsion
    :param increment: The grid increment, in radians
    :return: The coordinate grid
    """
    # [0, 30, 60, 90, ...] << in degrees
    vals = numpy.arange(0, span(tor), increment)
    # Start from the equilibrium value
    vals += zmat.value(zma, name(tor))
    return tuple(map(float, vals))


# # Transformations
def with_geometry_indices(tor: Torsion, zc_: ZmatConv) -> Torsion:
    """Given a z-matrix torsion and a z-matrix conversion, return a torsion with
    geometry indices.

    That is, the axis and group indices will skip dummy atoms

    :param tor: A torsion data structure
    :param zc_: A z-matrix conversion data structure
    :return: A torsion data structure using geometry indices
    """
    return from_data(
        name_=name(tor),
        coo=coordinate(tor, key_typ="geom", zc_=zc_),
        grps=groups(tor, key_typ="geom", zc_=zc_),
        symm=symmetry(tor),
        ngrps=neighbor_groups(tor, key_typ="geom", zc_=zc_),
    )


def update_against_zmatrix(tor: Torsion, zma: Any) -> Torsion:
    """Update a torsion object from a z-matrix.

    Temporarily needed to make sure the torsion object contains sufficient information

    :param tor: A torsion
    :param zma: A z-matrix to update against
    :return: The updated torsion
    """
    coo = list(reversed(zmat.coordinate(zma, name(tor))))
    # In case the coordinate was re-ordered, make sure the groups match
    grps = tuple(sorted(groups(tor), key=lambda g: coo[-1] in g))
    # Re-evaluate the neighbor groups, in case they are missing
    nkeys_dct = zmat.neighbor_keys(zma)
    ngrps = tuple(sorted(nkeys_dct[k]) for k in coo[1:3])
    return from_data(
        name_=name(tor), coo=coo, grps=grps, symm=symmetry(tor), ngrps=ngrps
    )


# Torsion list functions
def torsions_string(tor_lst: list[Torsion], one_indexed: bool = True) -> str:
    """Write a list of torsions to a string.

    :param tor_lst: A list of torsions
    :param one_indexed: Is this a one-indexed string? defaults to True
    :returns: A string representations of the torsions, as a flattened list
    """
    tor_yml_dct = torsions_yaml_data(tor_lst, one_indexed=one_indexed)
    tor_str = yaml.dump(tor_yml_dct, sort_keys=False)
    return tor_str


def torsions_from_string(tor_str: str, one_indexed: bool = True) -> list[Torsion]:
    """Write a list of torsions to a string.

    :param tor_str: A string representations of the torsions, as a flattened list
    :param one_indexed: Is this a one-indexed string? defaults to True
    :returns: A list of torsions
    """
    tor_yml_dct = yaml.load(tor_str, Loader=yaml.FullLoader)
    tor_lst = torsions_from_yaml_data(tor_yml_dct, one_indexed=one_indexed)
    return tor_lst


def torsions_yaml_data(tor_lst: list[Torsion], one_indexed: bool = True) -> dict:
    """Write a list of torsions to a yaml-formatted dictionary.

    :param tor_lst: A list of torsions
    :param one_indexed: Is this a one-indexed string? defaults to True
    :returns: A string representations of the torsions, as a flattened list
    """
    shift = 1 if one_indexed else 0

    tor_yml_dct = {}
    for tor in tor_lst:
        coo = [k + shift for k in coordinate(tor)]
        axis1, axis2 = (k + shift for k in axis(tor))
        coo_str = "-".join("*" if k is None else str(k) for k in coo)
        yml_dct = {
            "axis1": axis1,
            "axis2": axis2,
            "symmetry": symmetry(tor),
            "coordinate": coo_str,
        }
        yml_dct = _write_groups_to_yaml_dict(
            yml_dct, groups(tor), shift=shift, key_prefix="group"
        )
        yml_dct = _write_groups_to_yaml_dict(
            yml_dct, neighbor_groups(tor), shift=shift, key_prefix="neighbor_group"
        )
        tor_yml_dct[name(tor)] = yml_dct
    return tor_yml_dct


def torsions_from_yaml_data(tor_yml_dct: dict, one_indexed: bool = True) -> dict:
    """Read a list of torsions out of a yaml-formatted torsion dictionary.

    :param tor_lst: A list of torsions
    :param one_indexed: Is this a one-indexed string? defaults to True
    :returns: A string representations of the torsions, as a flattened list
    """
    shift = -1 if one_indexed else 0

    tor_lst = []
    for name_, yml_dct in tor_yml_dct.items():
        ax_ = list(map(yml_dct.get, ["axis1", "axis2"]))
        ax_ = [k + shift for k in ax_]
        grps = _read_groups_from_yaml_dict(yml_dct, shift, "group")
        ngrps = _read_groups_from_yaml_dict(yml_dct, shift, "neighbor_group")
        coo = yml_dct.get("coordinate", None)
        if coo is None:
            coo = [None, *ax_, None]
        else:
            coo = [None if s == "*" else int(s) + shift for s in coo.split("-")]

        tor = from_data(
            name_=name_,
            coo=coo,
            grps=grps,
            symm=yml_dct["symmetry"],
            ngrps=ngrps,
        )

        tor_lst.append(tor)

    return tuple(tor_lst)


# helper functions
def _write_groups_to_yaml_dict(
    yml_dct: dict[object, object],
    grps: Groups | None,
    shift: int,
    key_prefix: str,
) -> dict[object, object]:
    """Write rotational groups to a YAML dictionary.

    :param yml_dct: YAML dictionary
    :param grps: Rotational groups
    :param shift: The shift to apply to the keys
    :param key_prefix: Prefix for key
    :return: Updated YAML dictionary
    """
    if grps is None:
        return yml_dct

    yml_dct = yml_dct.copy()
    grp_vals = tuple(_group_to_yaml_value(g, shift) for g in grps)
    yml_dct.update({f"{key_prefix}{i+1}": v for i, v in enumerate(grp_vals)})
    return yml_dct


def _read_groups_from_yaml_dict(
    yml_dct: dict[object, object],
    shift: int,
    key_prefix: str,
) -> Groups | None:
    """Read rotational groups from a YAML dictionary.

    :param yml_dct: YAML dictionary
    :param shift: The shift to apply to the keys
    :param key_prefix: Prefix for key
    :return: Updated YAML dictionary
    """
    grp_vals = [yml_dct.get(f"{key_prefix}{i+1}") for i in range(2)]
    any_none = any(v is None for v in grp_vals)
    all_none = all(v is None for v in grp_vals)
    assert all_none or not any_none, f"In consistent group values: {yml_dct}"
    if any_none:
        return None
    return tuple(_group_from_yaml_value(v, shift) for v in grp_vals)


def _group_to_yaml_value(grp: Group, shift: int) -> str | int:
    """Get a YAML value to represent a group.

    :param grp: The rotational group
    :param shift: The shift to apply to the keys
    :return: The YAML value for the group
    """
    grp = [k + shift for k in grp]
    return "-".join(map(str, grp)) if len(grp) > 1 else grp[0]


def _group_from_yaml_value(grp_val: str | int, shift: int) -> Group:
    """Get a group from a YAML value.

    :param grp_str: The string representation of the rotational group
    :param shift: The shift to apply to the keys
    :return: The rotational group
    """
    grp = [grp_val] if isinstance(grp_val, int) else map(int, grp_val.split("-"))
    grp = [k + shift for k in grp]
    return grp
