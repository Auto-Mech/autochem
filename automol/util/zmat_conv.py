"""A data structure for tracking z-matrix conversion information

Format:

    zmat_conv = {
        0: (0, None, None),
        1: (1, None, None),
        2: (None, 1,    0), # < dummy off of atom 1, perpendicular to 0
        3: (2, None, None),
        4: (3, None, None),
        5: (None, 3,    2), # < dummy off of atom 3, perpendicular to 2
        6: (4, None, None),
        ...
        zkey: (gkey, lin_gkey, dir_gkey),
    }

    zkey is the z-matrix key for an atom
    gkey is the geometry key for an atom
        ('None' for dummy atoms)
    lin_gkey is the geometry key for the linear atom requiring a dummy atom
        ('None' for real atoms)
    dir_gkey is the geometry key for a linear atom neighbor, giving the linear direction
        ('None' for real atoms)

Note to self: The "direction" key should maybe be taken out. The one place where I was
using it, I don't think it serves a function anymore, and it gave the wrong choice for
the case that I was looking at (an Sn2 reaction).

Dummy atoms will be placed over the linear atom, perpendicular to the linear direction
as specified by the direction atom.

Together, the linear and direction keys will be referred to as the dummy atom's "source
keys", since they specify where the dummy atom is to be located.

Using the `relabel_dict` and `insert_dict` functions below, a z-matrix transformation
can be replicated in two steps:

    1. Relabel real atoms using the `relabel_dict`. Simply map geometry keys onto
    z-matrix keys:

        zkeys = [rel_dct[gkey] for gkey in gkeys]

    2. Insert dummy atoms using the `insert_dict`. The keys of this dictionary are the
    z-matrix keys of the dummy atoms, and the values are the z-matrix keys of their
    parent atoms.

        for dummy_zkey, parent_zkey in ins_dct.items():
            <Insert an atom connected to `parent_zkey` with key `dummy_zkey`>
"""
from typing import Dict, List, Optional, Tuple

from ._util import translate

Key = Optional[int]
ZmatConv = Dict[int, Tuple[Key, Key, Key]]


# Constructors
def from_zmat_data(zcount: int, src_zkeys_dct: Dict[int, Tuple[int, int]]) -> ZmatConv:
    """Generate a z-matrix conversion data structure from z-matrix information.

    :param zcount: The number of z-matrix atoms, including dummy atoms
    :param src_zkeys_dct: The source atoms for each dummy atom, by dummy key; the source
        atoms are (1.) the linear atom it connects to, and (2.) an atom specifying the
        linear direction, with respect to which it has a 90 degree angle
    :return: The z-matrix conversion.
    """
    # 1. Generate the mapping of final keys onto original keys
    real_zkeys = [k for k in range(zcount) if k not in src_zkeys_dct]
    gkey_dct = dict(map(reversed, enumerate(real_zkeys)))

    # 2. Apply the mapping to get the source atoms as geometry keys
    src_gkeys_dct = {
        k: tuple(map(gkey_dct.__getitem__, v)) for k, v in src_zkeys_dct.items()
    }

    # 3. Generate the z-matrix conversion
    zc_ = {zk: (gk, None, None) for zk, gk in gkey_dct.items()}
    zc_.update({zk: (None, *sgks) for zk, sgks in src_gkeys_dct.items()})

    # 4. Sort it for convenience
    zc_ = dict(sorted(zc_.items()))
    return zc_


def from_geom_data(gcount: List[int], lin_gkeys: List[int]) -> ZmatConv:
    """Generate a z-matrix conversion data structure from geometry information.

    The direction atoms are left unspecified (None)

    :param gcount: The number of geometry atoms, excluding dummy atoms
    :param lin_gkeys: The linear atom keys that require dummy atoms
    :return: The z-matrix conversion
    """
    assert all(
        0 <= k <= gcount for k in lin_gkeys
    ), f"Linear atom keys {lin_gkeys} are out of range [0, {gcount}]"

    # 1. Add real atoms to the list of values
    zc_vals = [(k, None, None) for k in range(gcount)]

    # 2. Add dummy atoms to the list of values
    zc_vals.extend((None, k, None) for k in lin_gkeys)

    # 3. Generate the z-matrix conversion
    zc_ = dict(enumerate(zc_vals))
    return zc_


# Getters
def count(zc_: ZmatConv, typ: str = "zmat") -> int:
    """Count the number of atoms.

    :param zc_: A z-matrix conversion data structure
    :param typ: The type of atoms to count ('geom' or 'zmat'); defaults to 'zmat'
    :return: The number of atoms
    """
    assert typ in ("geom", "zmat")

    return len(real_keys(zc_)) if typ == "geom" else len(zc_)


def real_keys(zc_: ZmatConv, typ: str = "zmat") -> List[int]:
    """Get the real atom keys from a z-matrix conversion

    :param zc_: A z-matrix conversion data structure
    :param typ: The type of keys to use ('geom' or 'zmat'); defaults to 'zmat'
    :return: The real atom keys
    """
    assert typ in ("geom", "zmat")

    keys = [
        gk if typ == "geom" else zk for zk, (gk, _, _) in zc_.items() if gk is not None
    ]
    return tuple(sorted(keys))


def dummy_keys(zc_: ZmatConv) -> List[int]:
    """Get the dummy atom keys from a z-matrix conversion

    :param zc_: A z-matrix conversion data structure
    :return: The dummy atom keys
    """
    keys = [zk for zk, (gk, _, _) in zc_.items() if gk is None]
    return tuple(sorted(keys))


def linear_atom_keys(zc_: ZmatConv, typ: str = "zmat") -> List[int]:
    """Get the list of original or final parent atom keyss from a z-matrix conversion.

    :param zc_: A z-matrix conversion data structure
    :param typ: The type of keys to use ('geom' or 'zmat'); defaults to 'zmat'
    :return: The list of parent atoms, using the original keys
    """
    assert typ in ("geom", "zmat")

    lin_keys = [lgk for _, (_, lgk, _) in zc_.items() if lgk is not None]
    if typ == "zmat":
        rel_dct = relabel_dict(zc_)
        lin_keys = tuple(map(rel_dct.__getitem__, lin_keys))
    return tuple(sorted(lin_keys))


def dummy_source_keys(zc_: ZmatConv, typ: str = "zmat") -> list[tuple[int, int]]:
    """Get the list of original or final parent atom keys from a z-matrix conversion.

    :param zc_: A z-matrix conversion data structure
    :param typ: The type of keys to use ('geom' or 'zmat'); defaults to 'zmat'
    :return: The list of parent atoms, using the original keys
    """
    assert typ in ("geom", "zmat")

    src_keys = [(lgk, dgk) for _, (_, lgk, dgk) in zc_.items() if lgk is not None]
    if typ == "zmat":
        rel_dct = relabel_dict(zc_)
        src_keys = tuple(tuple(map(rel_dct.__getitem__, ks)) for ks in src_keys)
    return tuple(sorted(src_keys))


def relabel_dict(zc_: ZmatConv, typ: str = "geom") -> Dict[int, int]:
    """Describes the relabeling of real (non-dummy) atoms.

    :param zc_: A z-matrix conversion data structure
    :param typ: The type of keys to be relabeled ('geom' or 'zmat'); defaults to 'geom'
    :return: A dictionary to relabel keys for dummy insertion.
    """  # noqa: D401
    assert typ in ("geom", "zmat")

    if typ == "geom":
        rel_dct = {gk: zk for zk, (gk, _, _) in zc_.items() if gk is not None}
    else:
        rel_dct = {zk: gk for zk, (gk, _, _) in zc_.items() if gk is not None}

    return rel_dct


def insert_dict(zc_: ZmatConv) -> Dict[int, int]:
    """Describes the insertion of dummies, mapping their keys upon insertion to
    their parent atom keys (*after relabeling*).

    :param zc_: A z-matrix conversion data structure
    :return: A dictionary to insert dummies after relabeling using `relabel_dict`.
    """  # noqa: D401
    rel_dct = relabel_dict(zc_)
    ins_dct = {zk: rel_dct[lgk] for zk, (gk, lgk, _) in zc_.items() if gk is None}
    return ins_dct


# Transformations
def without_direction_keys(zc_: ZmatConv) -> ZmatConv:
    """Remove direction keys from the z-matrix conversion.

    :param zc_: A z-matrix conversion data structure
    :return: A new z-matrix conversion data structure
    """
    zc_ = {zk: (gk, lgk, None) for zk, (gk, lgk, _) in zc_.items()}
    return zc_


def relabel(zc_: ZmatConv, key_dct: Dict[int, int], typ: str = "zmat") -> ZmatConv:
    """Relabel the final (or original) keys of a z-matrix conversion.

    :param zc_: A z-matrix conversion data structure
    :param key_dct: The relabeling
    :param typ: The type of keys to use ('geom' or 'zmat'); defaults to 'zmat'
    :return: The z-matrix conversion
    """
    assert typ in ("geom", "zmat")

    if typ == "geom":
        key_dct[None] = None
        zc_ = {zk: tuple(map(key_dct.__getitem__, gks)) for zk, gks in zc_.items()}
    else:
        zc_ = {key_dct[zk]: gks for zk, gks in zc_.items()}
    return zc_


def subset(
    zc_: ZmatConv, keys: List[int], typ: str = "zmat", dir_: Optional[bool] = None
) -> ZmatConv:
    """Extract the z-matrix conversion for a subset of keys

    :param zc_: A z-matrix conversion data structure
    :param keys: A subset of the keys in the conversion
    :param typ: The type of keys to use ('geom' or 'zmat'); defaults to 'zmat'
    :param dir_: Keep the direction keys? defaults to None, which keeps them for
        z-matrix subsets but not for geometry subsets
    :return: The subset z-matrix conversion
    """
    assert typ in ("geom", "zmat")

    if dir_ is False:
        zc_ = without_direction_keys(zc_)

    if typ == "geom":
        # If `keep_dir` is unspecified, drop direction keys
        zc_ = without_direction_keys(zc_) if dir_ is None else zc_

        # Grab the subset matching the geometry keys
        keys = list(keys) + [None]
        zc_ = {zk: gks for zk, gks in zc_.items() if all(k in keys for k in gks)}
    else:
        # Grab the subset matching the z-matrix keys
        zc_ = {zk: gks for zk, gks in zc_.items() if zk in keys}

    return zc_


# Comparisons
def isomorphism(
    zc1: ZmatConv, zc2: ZmatConv, sub: bool = False, dir_: bool = False
) -> Dict[int, int]:
    """Get the isomorphism of one ZmatConv onto another

    Assuming they share the same geometry keys, this backs out the graph isomorphism of
    one set of z-matrix keys onto another.

    :param zc1: A z-matrix conversion
    :type zc1: ZmatConv
    :param zc2: A relabeled z-matrix conversion
    :type zc2: ZmatConv
    :param sub: Find a subset isomorphism instead? defaults to False
    :type sub: bool, optional
    :param dir_: Include direction keys in the comparison? defaults to False
    :type dir: bool, optional
    :returns: A mapping of `zc1` onto `zc2`, if it exists
    """
    if not dir_:
        zc1 = without_direction_keys(zc1)
        zc2 = without_direction_keys(zc2)

    keys1 = list(zc1.keys())
    vals1 = list(zc1.values())
    keys2 = list(zc2.keys())
    vals2 = list(zc2.values())

    if sub:
        if not (len(zc1) >= len(zc2) and set(vals1) >= set(vals2)):
            return None
    else:
        if not (len(zc1) == len(zc2) and set(vals1) == set(vals2)):
            return None

    key_dct = {keys1[vals1.index(v2)]: k2 for k2, v2 in zip(keys2, vals2)}
    return key_dct


def subset_isomorphism(
    zc1: ZmatConv, zc2: ZmatConv, rev: bool = False, dir_: bool = False
) -> Dict[int, int]:
    """Get the isomorphism of a subset of one ZmatConv onto another

    Assuming they share the same geometry keys, this backs out the graph isomorphism of
    a subset of the z-matrix keys in `zc1` onto the z-matrix keys in `zc2`.

    :param zc1: A z-matrix conversion
    :type zc1: ZmatConv
    :param zc2: A relabeled z-matrix conversion
    :type zc2: ZmatConv
    :param rev: Get the reverse isomorphism, of `zc2` onto `zc1`? defaults to False
    :type rev: bool, optional
    :param dir_: Include direction keys in the comparison? defaults to False
    :type dir: bool, optional
    :returns: A mapping of `zc1` onto `zc2`, if it exists
    """
    iso_dct = isomorphism(zc1, zc2, sub=True, dir_=dir_)
    return dict(map(reversed, iso_dct.items())) if rev else iso_dct


# Utility functions
def relabel_zmatrix_key_sequence(
    zc_: ZmatConv, zkeys: List[int], dummy: bool = False
) -> List[int]:
    """Convert a (potentially nested) list of z-matrix keys into geometry keys

    Relabels real atoms and drops dummy atoms

    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :param zkeys: A list of z-matrix keys, nested to a consistent depth
    :type zkeys: List[int]
    :param dummy: Keep the dummy atoms as `None` values? defaults to False
    :type dummy: bool, optional
    :return: The corresponding keys of the geometry, minus dummy atoms
    :rtype: List[int]
    """
    rel_dct = relabel_dict(zc_, typ="zmat")
    return translate(zkeys, rel_dct, drop=not dummy)


# Conversions
def yaml_data(zc_: ZmatConv) -> list:
    """A yaml-friendly format for the z-matrix conversion

    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :return: A yaml-formatted z-matrix conversion
    :rtype: list
    """
    zc_yml = [[k, *v] for k, v in sorted(zc_.items())]
    return zc_yml


def from_yaml_data(zc_yml: list) -> ZmatConv:
    """Put a yaml-formatted z-matrix conversion back into standard format

    :param zc_yml: A yaml-formatted z-matrix conversion
    :type zc_yml: list
    :returns: A z-matrix conversion data structure
    :rtype: ZmatConv
    """
    keys = [row[0] for row in zc_yml]
    vals = [tuple(row[1:]) for row in zc_yml]
    zc_ = dict(zip(keys, vals))
    return zc_
