"""Functions for keeping track of dummy atoms using a "dummy transformation table"

Format:

    dummy_trans_table = {
        0: (0, None),
        1: (1, None),  # < parent atom 1
        2: (None, 1),  # < dummy atom off of parent atom 1
        3: (2, None),
        4: (3, None),  # < parent atom 2
        5: (None, 3),  # < dummy atom off of parent atom 2
        6: (4, None),
        ...
        k: (k0, k0p),
    }

    k is the new key for the atom (dummy or parent)
    k0 is the old key for the atom (`None` for dummy atoms)
    k0p is the old key of the parent atom (`None` for non-dummy atoms)

Using the `relabel_dict` and `insert_dict` functions below, dummy atoms can be inserted
using a two-step process:

    1. Relabel parent atoms using the result of `relabel_dict`. Simply map old keys onto
    new keys: new_key = rel_dct[old_key].
    2. Insert dummy atoms using the result of `insert_dict`. The keys are the new dummy
    atom keys, the values are the *new* parent atom keys (after relabeling):
        for dummy_key, parent_key in ins_dct.items():
            <Insert an atom connected to `parent_key` with key `dummy_key`>
"""
from typing import Dict, List, Tuple, Union

Key = Union[int, None]
DummyTransTableRow = Tuple[Key, Key, Key]
DummyTransTable = List[DummyTransTableRow]


def from_parent_atoms_list(
    nk0s: int, k0ps: List[int], insert: bool = False
) -> DummyTransTable:
    """Generate a DummyTransTable from a list of parent atoms

    :param nk0s: The original number of atoms
    :type nk0s: int
    :param k0ps: The parent atoms which will receive dummy atoms from the original list
    :type k0ps: List[int]
    :param insert: Insert the dummy atoms immediately after their parents? If not, they
        will be appended to the end; defaults to False
    :type insert: bool, optional
    :return: The dummy transformation table
    :rtype: DummyTransTable
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


def to_relabel_dict(dtt: DummyTransTable, rev: bool = False) -> Dict[int, int]:
    """Describes the relabeling of non-dummy atoms

    :param dtt: The dummy table
    :type dtt: DummyTransTable
    :param rev: Get the reverse mapping?, defaults to False
    :type rev: bool, optional
    :return: A dictionary to relabel keys for dummy atom insertion
    :rtype: Dict[int, int]
    """
    rel_dct = {k0: k for k, (k0, _) in dtt.items() if k0 is not None}
    if rev:
        rel_dct = dict(map(reversed, rel_dct.items()))

    return rel_dct


def to_insert_dict(dtt: DummyTransTable) -> Dict[int, int]:
    """Describes the insertion of dummy atoms, mapping their keys upon insertion to
    their parent atom keys (*after relabeling*)

    :param dtt: The dummy table
    :type dtt: DummyTransTable
    :return: A dictionary to insert dummy atoms after relabeling using `relabel_dict`
    :rtype: Dict[int, int]
    """
    rel_dct = to_relabel_dict(dtt)
    ins_dct = {k: rel_dct[k0p] for k, (k0, k0p) in dtt.items() if k0 is None}
    return ins_dct


def parent_atom_list(dtt: DummyTransTable) -> List[int]:
    """Get the list of parent atoms from a dummy transformation table, using the
    original keys

    :param dtt: The dummy table
    :type dtt: DummyTransTable
    :return: The list of parent atoms, using the original keys
    :rtype: List[int]
    """
    k0ps = tuple(k0p for _, (_, k0p) in dtt.items() if k0p is not None)
    return k0ps


def relabel(dtt: DummyTransTable, key_dct: Dict[int, int]) -> DummyTransTable:
    """Relabel the output keys of a dummy transformation

    :param dtt: The dummy table
    :type dtt: DummyTransTable
    :param key_dct: The relabeling
    :type key_dct: Dict[int, int]
    :return: The list of parent atoms, using the original keys
    :rtype: List[int]
    """
    dtt = {key_dct[k]: (k0, k0p) for k, (k0, k0p) in dtt.items()}
    return dtt


def isomorphism(dtt1: DummyTransTable, dtt2: DummyTransTable):
    """Get the isomorphism of on DummyTransTable onto another

    Assumes they both share the same original keys

    :param dtt1: A dummy transformation table
    :type dtt1: DummyTransTable
    :param dtt2: A relabeled dummy transformation table
    :type dtt2: DummyTransTable
    """
    assert (
        len(dtt1) == len(dtt2) and set(dtt1.values()) == set(dtt2.values())
    ), f"These DummyTransTables don't have the same original keys:\n{dtt1}\n{dtt2}"

    keys1 = list(dtt1.keys())
    vals1 = list(dtt1.values())
    keys2 = list(dtt2.keys())
    vals2 = list(dtt2.values())

    key_dct = {k1: keys2[vals2.index(v1)] for k1, v1 in zip(keys1, vals1)}
    return key_dct
