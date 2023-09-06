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


def from_original_parent_atom_keys(
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


def true_atom_keys(dc_: DummyConv, original: bool = False) -> List[int]:
    """Get the list of true atoms from a dummy conversion

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :param original: Return the original, instead of the final keys?
    :type original: bool
    :return: The list of parent atoms, using the original keys
    :rtype: List[int]
    """
    par_keys = [k0 if original else k for k, (k0, _) in dc_.items() if k0 is not None]
    return tuple(sorted(par_keys))


def dummy_atom_keys(dc_: DummyConv) -> List[int]:
    """Get the list of dummy atoms from a dummy conversion, using the
    final keys

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :return: The list of parent atoms, using the original keys
    :rtype: List[int]
    """
    par_keys = [k for k, (k0, _) in dc_.items() if k0 is None]
    return tuple(sorted(par_keys))


def parent_atom_keys(dc_: DummyConv, original: bool = False) -> List[int]:
    """Get the list of parent atoms from a dummy conversion, using the
    original keys

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


def relabel(dc_: DummyConv, key_dct: Dict[int, int]) -> DummyConv:
    """Relabel the output keys of a dummy conversion

    :param dc_: A dummy conversion data structure
    :type dc_: DummyConv
    :param key_dct: The relabeling
    :type key_dct: Dict[int, int]
    :return: The dummy conversion
    :rtype: DummyConv
    """
    dc_ = {key_dct[k]: (k0, k0p) for k, (k0, k0p) in dc_.items()}
    return dc_


def isomorphism(dc1: DummyConv, dc2: DummyConv):
    """Get the isomorphism of on DummyConv onto another

    Assumes they both share the same original keys

    :param dc1: A dummy conversion
    :type dc1: DummyConv
    :param dc2: A relabeled dummy conversion
    :type dc2: DummyConv
    """
    assert len(dc1) == len(dc2) and set(dc1.values()) == set(
        dc2.values()
    ), f"These DummyConvs don't have the same original keys:\n{dc1}\n{dc2}"

    keys1 = list(dc1.keys())
    vals1 = list(dc1.values())
    keys2 = list(dc2.keys())
    vals2 = list(dc2.values())

    key_dct = {k1: keys2[vals2.index(v1)] for k1, v1 in zip(keys1, vals1)}
    return key_dct
