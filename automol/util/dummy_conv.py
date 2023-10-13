"""A data structure for keeping track of conversions involving dummy atoms

Format:

    dummy_conv = {
        0: (0, None),
        1: (1, None),  # < parent atom 1
        2: (None, 1),  # < dummy off of parent atom 1
        3: (2, None),
        4: (3, None),  # < parent atom 2
        5: (None, 3),  # < dummy off of parent atom 2
        6: (4, None),
        ...
        k: (k0, k0p),
    }

    k is the new key for the atom (dummy or parent)
    k0 is the old key for the atom (`None` for dummies)
    k0p is the old key of the parent atom (`None` for non-dummies)

Using the `relabel_dict` and `insert_dict` functions below, dummies can be inserted
using a two-step process:

    1. Relabel parent atoms using the result of `relabel_dict`. Simply map old keys onto
    new keys: new_key = rel_dct[old_key].
    2. Insert dummies using the result of `insert_dict`. The keys are the new dummy
    atom keys, the values are the *new* parent atom keys (after relabeling):
        for dummy_key, parent_key in ins_dct.items():
            <Insert an atom connected to `parent_key` with key `dummy_key`>
"""
from typing import Dict, List, Tuple, Union

Key = Union[int, None]
DummyConv = Dict[int, Tuple[Key, Key]]


def from_original_dummy_parent_keys(
    nk0s: int, k0ps: List[int], insert: bool = False
) -> DummyConv:
    """Generate a dummy conversion data structure from a list of original parent keys

    :param nk0s: The original number of atoms
    :type nk0s: int
    :param k0ps: The parent atoms which will receive dummies from the original list
    :type k0ps: List[int]
    :param insert: Insert the dummies immediately after their parents? If not, they
        will be appended to the end; defaults to False
    :type insert: bool, optional
    :return: The dummy conversion
    :rtype: DummyConv
    """
    k0s = list(range(nk0s))
    vals = k0s.copy()

    for k0p in k0ps:
        if insert:
            vals.insert(vals.index(k0p) + 1, (None, k0p))
        else:
            vals.append((None, k0p))

    vals = [(val, None) if val in k0s else val for val in vals]

    return dict(enumerate(vals))


def from_insert_dict(nks: int, kp_dct: Dict[int, int]) -> DummyConv:
    """Generate a dummy conversion data structure from an insertion dictionary

    (Final parent atom keys, by dummy atom key)

    :param nks: The final number of atoms
    :type nks: int
    :param kp_dct: Final parent keys, by final dummy key
    :type kp_dct: Dict[int, int]
    :return: The dummy conversion
    :rtype: DummyConv
    """
    # 1. Generate the mapping of final keys onto original keys
    k0_dct = dict(map(reversed, enumerate(k for k in range(nks) if k not in kp_dct)))
    # 2. Generate the dummy conversion
    dc_ = {
        k: (None, k0_dct[kp_dct[k]]) if k in kp_dct else (k0_dct[k], None)
        for k in range(nks)
    }
    return dc_


def relabel_dict(dc_: DummyConv, rev: bool = False) -> Dict[int, int]:
    """Describes the relabeling of non-dummies

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :param rev: Get the reverse mapping?, defaults to False
    :type rev: bool, optional
    :return: A dictionary to relabel keys for dummy insertion
    :rtype: Dict[int, int]
    """
    rel_dct = {k0: k for k, (k0, _) in dc_.items() if k0 is not None}
    if rev:
        rel_dct = dict(map(reversed, rel_dct.items()))

    return rel_dct


def insert_dict(dc_: DummyConv) -> Dict[int, int]:
    """Describes the insertion of dummies, mapping their keys upon insertion to
    their parent atom keys (*after relabeling*)

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :return: A dictionary to insert dummies after relabeling using `relabel_dict`
    :rtype: Dict[int, int]
    """
    rel_dct = relabel_dict(dc_)
    ins_dct = {k: rel_dct[k0p] for k, (k0, k0p) in dc_.items() if k0 is None}
    return ins_dct


def count(dc_: DummyConv, original: bool = False) -> int:
    """Count the number of atoms

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :param original: Count the original, instead of the final atoms?
    :type original: bool
    :return: The number of atoms
    :rtype: int
    """
    return len(real_keys(dc_)) if original else len(dc_)


def real_keys(dc_: DummyConv, original: bool = False) -> List[int]:
    """Get the list of real atoms from a dummy conversion

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :param original: Return the original, instead of the final keys?
    :type original: bool
    :return: The list of parent atoms, using the original keys
    :rtype: List[int]
    """
    par_keys = [k0 if original else k for k, (k0, _) in dc_.items() if k0 is not None]
    return tuple(sorted(par_keys))


def dummy_keys(dc_: DummyConv) -> List[int]:
    """Get the list of dummy atoms from a dummy conversion, using the
    final keys

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :return: The list of parent atoms, using the original keys
    :rtype: List[int]
    """
    dum_keys = [k for k, (k0, _) in dc_.items() if k0 is None]
    return tuple(sorted(dum_keys))


def dummy_parent_keys(dc_: DummyConv, original: bool = False) -> List[int]:
    """Get the list of original or final parent atom keyss from a dummy conversion

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :param original: Return the original, instead of the final keys?
    :type original: bool
    :return: The list of parent atoms, using the original keys
    :rtype: List[int]
    """
    par_keys = [k0p for _, (_, k0p) in dc_.items() if k0p is not None]
    if not original:
        rel_dct = relabel_dict(dc_)
        par_keys = tuple(map(rel_dct.__getitem__, par_keys))
    return tuple(sorted(par_keys))


def shift(dc_: DummyConv, start: int, original: bool = False) -> DummyConv:
    """Shift the final (or original) keys of a dummy conversion to start from a
    different value

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :param start: This will be the new starting index
    :type start: int
    :param original: Shift the original, instead of the final keys? defaults to False
    :type original: bool, optional
    :return: The dummy conversion
    :rtype: DummyConv
    """
    keys0 = sorted(real_keys(dc_, original=True) if original else dc_.keys())
    key_dct = {k: k + start for k in keys0}
    return relabel(dc_, key_dct, original=original)


def relabel(
    dc_: DummyConv, key_dct: Dict[int, int], original: bool = False
) -> DummyConv:
    """Relabel the final (or original) keys of a dummy conversion

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :param key_dct: The relabeling
    :type key_dct: Dict[int, int]
    :param original: Relabel the original, instead of the final keys? default False
    :type original: bool, optional
    :return: The dummy conversion
    :rtype: DummyConv
    """
    if original:
        key_dct[None] = None
        dc_ = {k: (key_dct[k0], key_dct[k0p]) for k, (k0, k0p) in dc_.items()}
    else:
        dc_ = {key_dct[k]: (k0, k0p) for k, (k0, k0p) in dc_.items()}
    return dc_


def isomorphism(
    dc1: DummyConv, dc2: DummyConv, subgraph: bool = False
) -> Dict[int, int]:
    """Get the isomorphism of one DummyConv onto another

    Assuming they share the same keys, this backs out the isomorphism of the graph
    corresponding to (the original keys of) `dc1` onto the graph corresponding to (the
    original keys of) `dc2`.

    :param dc1: A dummy conversion
    :type dc1: DummyConv
    :param dc2: A relabeled dummy conversion
    :type dc2: DummyConv
    :param subgraph: Find a subgraph isomorphism instead? defaults to False
    :type subgraph: bool, optional
    :returns: A mapping of `dc1` onto `dc2`, if it exists
    """
    keys1 = list(dc1.keys())
    vals1 = list(dc1.values())
    keys2 = list(dc2.keys())
    vals2 = list(dc2.values())

    if subgraph:
        if not (len(dc1) >= len(dc2) and set(vals1) >= set(vals2)):
            return None
    else:
        if not (len(dc1) == len(dc2) and set(vals1) == set(vals2)):
            return None

    key_dct = {keys1[vals1.index(v2)]: k2 for k2, v2 in zip(keys2, vals2)}
    return key_dct


def subgraph_isomorphism(
    dc1: DummyConv, dc2: DummyConv, reverse: bool = False
) -> Dict[int, int]:
    """Get the isomorphism of a subset of one DummyConv onto another

    Assuming they share the same keys, this backs out the isomorphism of a subgraph of
    the graph corresponding to (the original keys of) `dc1` onto the graph corresponding
    to (the original keys of) `dc2`.

    :param dc1: A dummy conversion
    :type dc1: DummyConv
    :param dc2: A relabeled dummy conversion
    :type dc2: DummyConv
    :param reverse: Get the reverse isomorphism, of `dc2` onto `dc1`? defaults to False
    :type reverse: bool, optional
    :returns: A mapping of `dc1` onto `dc2`, if it exists
    """
    iso_dct = isomorphism(dc1, dc2, subgraph=True)
    return dict(map(reversed, iso_dct.items())) if reverse else iso_dct


def yaml_data(dc_: DummyConv) -> list:
    """A yaml-friendly format for the dummy conversion

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :return: A yaml-formatted dummy conversion
    :rtype: list
    """
    dc_yml = [[k, *v] for k, v in sorted(dc_.items())]
    return dc_yml


def from_yaml_data(dc_yml: list) -> DummyConv:
    """Put a yaml-formatted dummy conversion back into standard format

    :param dc_yml: A yaml-formatted dummy conversion
    :type dc_yml: list
    :returns: A dummy conversion data structure
    :rtype: DummyConv
    """
    keys = [row[0] for row in dc_yml]
    vals = [tuple(row[1:]) for row in dc_yml]
    dc_ = dict(zip(keys, vals))
    return dc_
