"""TS classification and other functions
"""

import functools
import itertools
from typing import Dict, Optional, Tuple

from ... import util
from ._00core import (
    AtomKey,
    BondKey,
    tetrahedral_atom_keys,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reactants_graph_without_stereo,
    ts_reverse,
    vinyl_radical_bond_candidates,
)
from ._02algo import reacting_rings_bond_keys


# reaction site classification
def atom_transfers(tsg) -> Dict[int, Tuple[int, int]]:
    """Get a dictionary describing atom transfers; keys are transferring atoms, values
    are donors and acceptors, respectively

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A list of triples containing the donor atom, the transferring atom, and
        the acceptor atom, respectively
    :rtype: Dict[int, Tuple[int, int]]
    """
    brk_bkeys = ts_breaking_bond_keys(tsg)
    frm_bkeys = ts_forming_bond_keys(tsg)

    tra_dct = {}
    for brk_bkey, frm_bkey in itertools.product(brk_bkeys, frm_bkeys):
        if brk_bkey & frm_bkey:
            (tra_key,) = brk_bkey & frm_bkey
            (don_key,) = brk_bkey - frm_bkey
            (acc_key,) = frm_bkey - brk_bkey
            tra_dct[tra_key] = (don_key, acc_key)

    return tra_dct


def substitutions(tsg) -> Dict[int, Tuple[int, int]]:
    """Get a dictionary describing substitution reaction sites

    Maps transferring atoms onto their leaving and entering atoms, respectively

    (Limited to substitutions at tetrahedral atoms)

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A mapping of transferring atoms onto leaving and entering atoms
    :rtype: Dict[int, Tuple[int, int]]
    """
    tra_dct = atom_transfers(tsg)
    tra_keys = set(tra_dct.keys())
    tet_keys = tetrahedral_atom_keys(tsg)
    subst_keys = tra_keys & tet_keys
    return util.dict_.by_key(tra_dct, subst_keys)


def vinyl_addition_candidates(
    tsg, min_ncount: int = 1
) -> Dict[BondKey, Tuple[AtomKey, AtomKey]]:
    """Get a dictionary describing vinyl addition reaction site *candidates*

    Maps vinyl bond keys onto vinyl radical atoms entering atoms, respectively

    (Only finds candidates, since resonance evaluation is expensive)

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param min_ncount: Minimal # neighbor keys for consideration
    :type min_ncount: int, optional
    :returns: A mapping of transferring atoms onto leaving and entering atoms
    """
    rcts_gra = ts_reactants_graph_without_stereo(tsg)

    frm_bkeys = ts_forming_bond_keys(tsg)
    vin_dct = vinyl_radical_bond_candidates(rcts_gra, min_ncount=min_ncount)

    vin_add_dct = {}
    for bkey, key in vin_dct.items():
        # Find a forming bond at the vinyl radical site
        frm_bkey = next((bk for bk in frm_bkeys if key in bk), None)
        if frm_bkey is not None:
            (ent_key,) = frm_bkey - {key}

            vin_add_dct[bkey] = (key, ent_key)
    return vin_add_dct


def eliminations(tsg) -> Dict[BondKey, Tuple[AtomKey, AtomKey, Optional[BondKey]]]:
    """Get a dictionary describing elimination reaction sites

    Maps bonds across which eliminations occur onto their leaving atoms, along with the
    forming bond key for each, if present

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A mapping of elimination bonds onto leaving atoms and forming bond keys
        (Leaving atoms are sorted in order of the elimination bond atoms)
    """
    brk_bkeys_pool = ts_breaking_bond_keys(tsg)
    frm_bkeys_pool = ts_forming_bond_keys(tsg)
    rng_bkeys_lst = reacting_rings_bond_keys(tsg)

    def is_elimination_bond_(brk_bkeys, frm_bkeys):
        """An elimination bond key is a ring key that intersects two breaking bonds and
        is not a forming bond
        """

        def is_elimination_bkey(bkey):
            return (
                all(bkey & bk for bk in brk_bkeys) and bkey not in frm_bkeys | brk_bkeys
            )

        return is_elimination_bkey

    def common_atom(bkey1, bkey2):
        return next(iter(bkey1 & bkey2))

    # 1. Check reacting rings
    elim_dct = {}
    for rng_bkeys in rng_bkeys_lst:
        brk_bkeys = brk_bkeys_pool & rng_bkeys
        frm_bkeys = frm_bkeys_pool & rng_bkeys

        # 2. Require two breaking bonds within the ring
        if len(brk_bkeys) == 2:
            # 3. Find the elimination bond
            is_elim = is_elimination_bond_(brk_bkeys, frm_bkeys)
            bkey = next(filter(is_elim, rng_bkeys), None)
            if bkey is not None:
                # a. Sort the breaking bonds in order of the elimination bond atom
                brk_bkey1, brk_bkey2 = sorted(
                    brk_bkeys, key=functools.partial(common_atom, bkey)
                )

                # b. Get the corresponding leaving atoms
                (lea_key1,) = brk_bkey1 - bkey
                (lea_key2,) = brk_bkey2 - bkey

                # c. Get the forming bond key, if there is one
                assert len(frm_bkeys) <= 1, "Unexpected multiple forming bonds:{tsg}"
                frm_bkey = next(iter(frm_bkeys), None)

                elim_dct[bkey] = (lea_key1, lea_key2, frm_bkey)

    return elim_dct


def insertions(tsg) -> Dict[int, Tuple[int, int]]:
    """Get a dictionary describing insertion reaction sites

    Maps bonds across which insertions occur onto their entering atoms, along with the
    breaking bond key for each, if present

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A mapping of insertion bonds onto entering atoms and breaking bond keys
        (Entering atoms are sorted in order of the insertion bond atoms)
    :rtype: Dict[int, Tuple[int, int]]
    """
    return eliminations(ts_reverse(tsg))


def ring_forming_scissions(tsg) -> Dict[int, Tuple[int, int]]:
    """Get a dictionary describing substitution reaction sites

    Maps transferring atoms onto their leaving and entering atoms, respectively

    (Limited to substitutions at tetrahedral atoms)

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :returns: A mapping of transferring atoms onto leaving and entering atoms
    :rtype: Dict[int, Tuple[int, int]]
    """
    tra_dct = atom_transfers(tsg)
    rngs_bkeys = reacting_rings_bond_keys(tsg)

    rsciss_dct = {}
    for tra_key, (brk_nkey, frm_nkey) in tra_dct.items():
        brk_bkey = frozenset({tra_key, brk_nkey})
        frm_bkey = frozenset({tra_key, frm_nkey})
        rng_bkeys = next(
            (r for r in rngs_bkeys if frm_bkey in r and brk_bkey not in r), None
        )
        if rng_bkeys is not None:
            rsciss_dct[tra_key] = (brk_nkey, frm_nkey)

    return rsciss_dct
