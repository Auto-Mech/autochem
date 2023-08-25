"""Functions for keeping track of dummy atoms using a "dummy table"

Format:

    dummy_table = [
        (0, 0, None),
        (1, 1, None),  # < parent atom 1
        (2, None, 1),  # < dummy atom off of parent atom 1
        (3, 2, None),
        (4, 3, None),  # < parent atom 2
        (5, None, 3),  # < dummy atom off of parent atom 2
        (6, 4, None),
        ...
        (k, k0, k0p),
    ]

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
from typing import Dict, List, NoneType, Tuple, Union

Key = Union[int, NoneType]
DummyTableRow = Tuple[Key, Key, Key]
DummyTable = List[DummyTableRow]


def relabel_dict(dt_: DummyTable, rev: bool = False) -> Dict[int, int]:
    """Describes the relabeling of non-dummy atoms

    :param dt_: The dummy table
    :type dt_: DummyTable
    :param rev: Get the reverse mapping?, defaults to False
    :type rev: bool, optional
    :return: A dictionary to relabel keys for dummy atom insertion
    :rtype: Dict[int, int]
    """
    rel_dct = {k0: k for k, k0, _ in dt_ if k0 is not None}
    if rev:
        rel_dct = dict(map(reversed, rel_dct.items()))

    return rel_dct


def insert_dict(dt_: DummyTable) -> Dict[int, int]:
    """Describes the insertion of dummy atoms, mapping their keys upon insertion to
    their parent atom keys (*after relabeling*)

    :param dt_: The dummy table
    :type dt_: DummyTable
    :return: A dictionary to insert dummy atoms after relabeling using `relabel_dict`
    :rtype: Dict[int, int]
    """
    rel_dct = relabel_dict(dt_)
    ins_dct = {k: rel_dct[k0p] for k, k0, k0p in dt_ if k0 is None}
    return ins_dct
