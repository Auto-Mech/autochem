"""A data structure for tracking z-matrix conversion information

Format:

    zmat_conv = {
        0: (0, None),
        1: (1, None),  # < parent atom 1
        2: (None, 1),  # < dummy off of parent atom 1
        3: (2, None),
        4: (3, None),  # < parent atom 2
        5: (None, 3),  # < dummy off of parent atom 2
        6: (4, None),
        ...
        zkey: (gkey, parent_gkey),
    }

    zkey is the z-matrix key for an atom
    gkey is the geometry key for an atom
        ('None' for dummy atoms)
    parent_gkey is the geometry key for the parent of a dummy atom
        ('None' for real atoms)

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
import numbers
from collections.abc import Sequence
from typing import Dict, List, Tuple, Union

Key = Union[int, None]
ZmatConv = Dict[int, Tuple[Key, Key]]


# Constructors
def from_geom_dummy_parent_keys(
    nreal: int, parent_gkeys: List[int], insert: bool = False
) -> ZmatConv:
    """Generate a z-matrix conversion data structure from a list of geometry parent keys

    :param nreal: The number of real atoms in the geometry
    :type nreal: int
    :param parent_gkeys: The geometry keys of parent atoms (atoms receiving dummy atoms)
    :type parent_gkeys: List[int]
    :param insert: Insert the dummies immediately after their parents? If not, they
        will be appended to the end; defaults to False
    :type insert: bool, optional
    :return: The z-matrix conversion
    :rtype: ZmatConv
    """
    gkeys = list(range(nreal))
    zc_vals = gkeys.copy()

    for parent_gkey in parent_gkeys:
        if insert:
            zc_vals.insert(zc_vals.index(parent_gkey) + 1, (None, parent_gkey))
        else:
            zc_vals.append((None, parent_gkey))

    zc_vals = [(val, None) if val in gkeys else val for val in zc_vals]

    return dict(enumerate(zc_vals))


def from_zmat_dummy_parent_dict(
    ntotal: int, parent_zkey_dct: Dict[int, int]
) -> ZmatConv:
    """Generate a z-matrix conversion data structure from a dictionary mapping dummy
    atom keys into their parent atom keys

    :param ntotal: The total number of atoms in the z-matrix, including dummy atoms
    :type ntotal: int
    :param parent_zkey_dct: Final parent keys, by final dummy key
    :type parent_zkey_dct: Dict[int, int]
    :return: The z-matrix conversion
    :rtype: ZmatConv
    """
    # 1. Generate the mapping of final keys onto original keys
    gkey_dct = dict(
        map(reversed, enumerate(k for k in range(ntotal) if k not in parent_zkey_dct))
    )
    # 2. Generate the z-matrix conversion
    zc_ = {
        zk: (None, gkey_dct[parent_zkey_dct[zk]])
        if zk in parent_zkey_dct
        else (gkey_dct[zk], None)
        for zk in range(ntotal)
    }
    return zc_


# Conversion functions
def relabel_dict(zc_: ZmatConv, rev: bool = False) -> Dict[int, int]:
    """Describes the relabeling of non-dummies

    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :param rev: Get the reverse mapping?, defaults to False
    :type rev: bool, optional
    :return: A dictionary to relabel keys for dummy insertion
    :rtype: Dict[int, int]
    """
    rel_dct = {gk: zk for zk, (gk, _) in zc_.items() if gk is not None}
    if rev:
        rel_dct = dict(map(reversed, rel_dct.items()))

    return rel_dct


def insert_dict(zc_: ZmatConv) -> Dict[int, int]:
    """Describes the insertion of dummies, mapping their keys upon insertion to
    their parent atom keys (*after relabeling*)

    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :return: A dictionary to insert dummies after relabeling using `relabel_dict`
    :rtype: Dict[int, int]
    """
    rel_dct = relabel_dict(zc_)
    ins_dct = {zk: rel_dct[pgk] for zk, (gk, pgk) in zc_.items() if gk is None}
    return ins_dct


def geometry_keys(zc_: ZmatConv, zkeys: List[int], dummy: bool = False) -> List[int]:
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
    rel_dct = relabel_dict(zc_, rev=True)

    def zkeys_to_gkeys_(zks):
        """Recursively convert a nested list of z-matrix keys to geometry keys"""

        assert isinstance(zks, Sequence), f"Cannot process non-sequence {zks}"

        # If the first element is an integer, assume we have a list of z-matrix keys and
        # conver them
        if isinstance(zks[0], numbers.Integral):
            if dummy:
                gks = tuple(rel_dct[k] if k in rel_dct else None for k in zks)
            else:
                gks = tuple(rel_dct[k] for k in zks if k in rel_dct)
        # If we have a sequence, recursively call this function
        elif isinstance(zks[0], Sequence):
            gks = tuple(map(zkeys_to_gkeys_, zks))

        return gks

    return zkeys_to_gkeys_(zkeys)


# Getters
def count(zc_: ZmatConv, typ: str = "zmat") -> int:
    """Count the number of atoms

    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :param typ: The type of atoms to count ('geom' or 'zmat'); defaults to 'zmat'
    :type typ: str, optional
    :return: The number of atoms
    :rtype: int
    """
    assert typ in ("geom", "zmat")

    return len(real_keys(zc_)) if typ == "geom" else len(zc_)


def real_keys(zc_: ZmatConv, typ: str = "zmat") -> List[int]:
    """Get the real atom keys from a z-matrix conversion

    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :param typ: The type of keys to use ('geom' or 'zmat'); defaults to 'zmat'
    :type typ: str, optional
    :return: The real atom keys
    :rtype: List[int]
    """
    assert typ in ("geom", "zmat")

    keys = [
        gk if typ == "geom" else zk for zk, (gk, _) in zc_.items() if gk is not None
    ]
    return tuple(sorted(keys))


def dummy_keys(zc_: ZmatConv) -> List[int]:
    """Get the dummy atom keys from a z-matrix conversion

    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :return: The dummy atom keys
    :rtype: List[int]
    """
    keys = [zk for zk, (gk, _) in zc_.items() if gk is None]
    return tuple(sorted(keys))


def dummy_parent_keys(zc_: ZmatConv, typ: str = "zmat") -> List[int]:
    """Get the list of original or final parent atom keyss from a z-matrix conversion

    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :param typ: The type of keys to use ('geom' or 'zmat'); defaults to 'zmat'
    :type typ: str, optional
    :return: The list of parent atoms, using the original keys
    :rtype: List[int]
    """
    assert typ in ("geom", "zmat")

    par_keys = [pgk for _, (_, pgk) in zc_.items() if pgk is not None]
    if typ == "zmat":
        rel_dct = relabel_dict(zc_)
        par_keys = tuple(map(rel_dct.__getitem__, par_keys))
    return tuple(sorted(par_keys))


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


# Other
def relabel(zc_: ZmatConv, key_dct: Dict[int, int], typ: str = "zmat") -> ZmatConv:
    """Relabel the final (or original) keys of a z-matrix conversion

    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :param key_dct: The relabeling
    :type key_dct: Dict[int, int]
    :param typ: The type of keys to use ('geom' or 'zmat'); defaults to 'zmat'
    :type typ: str, optional
    :return: The z-matrix conversion
    :rtype: ZmatConv
    """
    assert typ in ("geom", "zmat")

    if typ == "geom":
        key_dct[None] = None
        zc_ = {zk: (key_dct[gk], key_dct[pgk]) for zk, (gk, pgk) in zc_.items()}
    else:
        zc_ = {key_dct[zk]: (gk, pgk) for zk, (gk, pgk) in zc_.items()}
    return zc_


def subset(zc_: ZmatConv, keys: List[int], typ: str = "zmat") -> ZmatConv:
    """Extract the z-matrix conversion for a subset of keys

    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :param keys: A subset of the keys in the conversion
    :type keys: List[int]
    :param typ: , defaults to "zmat"
    :type typ: str, optional
    :return: The subset z-matrix conversion
    :rtype: ZmatConv
    """
    assert typ in ("geom", "zmat")

    if typ == "geom":
        zc_ = {
            zk: (gk, pgk) for zk, (gk, pgk) in zc_.items() if gk in keys or pgk in keys
        }
    else:
        zc_ = {zk: (gk, pgk) for zk, (gk, pgk) in zc_.items() if zk in keys}
    return zc_


def isomorphism(zc1: ZmatConv, zc2: ZmatConv, sub: bool = False) -> Dict[int, int]:
    """Get the isomorphism of one ZmatConv onto another

    Assuming they share the same geometry keys, this backs out the graph isomorphism of
    one set of z-matrix keys onto another.

    :param zc1: A z-matrix conversion
    :type zc1: ZmatConv
    :param zc2: A relabeled z-matrix conversion
    :type zc2: ZmatConv
    :param sub: Find a subset isomorphism instead? defaults to False
    :type sub: bool, optional
    :returns: A mapping of `zc1` onto `zc2`, if it exists
    """
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
    zc1: ZmatConv, zc2: ZmatConv, rev: bool = False
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
    :returns: A mapping of `zc1` onto `zc2`, if it exists
    """
    iso_dct = isomorphism(zc1, zc2, sub=True)
    return dict(map(reversed, iso_dct.items())) if rev else iso_dct
