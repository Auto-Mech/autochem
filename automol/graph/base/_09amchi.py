""" AutoMechanic Chemical Identifier (AMChI) generating functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools
import re
from typing import Dict, List

from ... import form, util
from ...amchi import base as amchi_base
from ._00core import (
    atom_implicit_hydrogens,
    atom_keys,
    atom_stereo_keys,
    atom_stereo_parities,
    atoms_neighbor_atom_keys,
    bond_stereo_keys,
    bond_stereo_parities,
    formula,
    relabel,
    terminal_atom_keys,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
)
from ._02algo import (
    connected_components,
    dfs_,
    dfs_children,
    dfs_missing_children,
    dfs_parents,
    rings_atom_keys,
)
from ._03kekule import (
    kekules_bond_orders_collated,
    vinyl_radical_atom_keys,
)
from ._08canon import canonical_amchi_graph_with_numbers


# AMChI functions
def amchi(gra, stereo: bool = True) -> str:
    """Generate an AMChI string with AMChI canonical numbers

    :param gra: A graph
    :type gra: automol graph data structure
    :param stereo: Include stereo, if present?, defaults to True
    :type stereo: bool, optional
    :returns: The AMChI string, along with the canonical numbers for each connected
    :rtype: str
    """
    chi, _ = amchi_with_numbers(gra, stereo=stereo)
    return chi


def amchi_with_numbers(gra, stereo: bool = True) -> (str, List[Dict[int, int]]):
    """Generate an AMChI string with canonical atom numbers

    :param gra: A graph
    :type gra: automol graph data structure
    :param stereo: Include stereo, if present?, defaults to True
    :type stereo: bool, optional
    :returns: The AMChI string, along with the canonical numbers for each connected
        component, as a list of dictionaries
    :rtype: (str, List[Dict[int, int]])
    """
    gras = connected_components(gra)
    chis, num_dcts = zip(*(_amchi_with_numbers(g, stereo=stereo) for g in gras))
    chi, srt = amchi_base.sorted_join(chis)
    num_dcts = tuple(num_dcts[i] for i in srt)
    return chi, num_dcts


def _amchi_with_numbers(gra, stereo=True):
    """single-component AMChI string from a connected graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereo in the AMChI string, if present?
    :type stereo: bool
    :returns: the AMChI string
    :rtype: str
    """
    # 1. Put the graph in canonical AMChI format
    gra, num_dct, is_can_dir, is_can_enant = canonical_amchi_graph_with_numbers(
        gra, stereo=stereo
    )

    # 3. Generate the appropriate layers
    #   a. Formula string
    fml_str = formula_string(gra)

    #   b. Main layers (c, h)
    main_lyr_dct = {"c": connection_layer(gra), "h": hydrogen_layer(gra)}

    #   b. Stereo layers (b, t, m, s)
    ste_lyr_dct = {"b": bond_stereo_layer(gra), "t": atom_stereo_layer(gra)}
    if is_can_enant is not None:
        ste_lyr_dct["m"] = "0" if is_can_enant else "1"
        ste_lyr_dct["s"] = "1"

    #   c. TS layers (k, f, r)
    ts_lyr_dct = {"k": bond_breaking_layer(gra), "f": bond_forming_layer(gra)}
    if is_can_dir is not None:
        ts_lyr_dct["r"] = "0" if is_can_dir else "1"

    # 4. Build the AMChI string
    chi = amchi_base.from_data(
        fml_lyr=fml_str,
        main_lyr_dct=main_lyr_dct,
        ste_lyr_dct=ste_lyr_dct,
        ts_lyr_dct=ts_lyr_dct,
    )

    return chi, num_dct


# # inchi checker
def inchi_is_bad(gra, ich):
    """Check if this is a bad InChI by comparing it to a graph with stereo

    The InChI is currently considered 'bad' if it:
    1. Is missing stereo
    2. Has mobile hydrogens

    :param gra: molecular graph, with stereo
    :type gra: automol graph data structure
    :param chi: ChI string
    :type chi: str
    :returns: True if the InChI is bad, otherwise False
    :rtype: bool
    """
    ste_atm_keys = atom_stereo_keys(gra)
    ste_bnd_keys = bond_stereo_keys(gra)

    is_bad = False

    # 1. Stereogenic nitrogen atoms
    nit_atm_keys = atom_keys(gra, symb="N")
    has_ste_nit_atm = bool(nit_atm_keys & ste_atm_keys)

    # 2. Stereogenic resonance bonds
    res_bnd_ords_dct = kekules_bond_orders_collated(gra)
    has_ste_res_bnd = any(1 in res_bnd_ords_dct[k] for k in ste_bnd_keys)

    # 3. Stereogenic vinyl radical bonds
    vin_atm_keys = vinyl_radical_atom_keys(gra)
    has_ste_vin_bnd = any(k & vin_atm_keys for k in ste_bnd_keys)

    is_bad |= (
        has_ste_nit_atm
        or has_ste_res_bnd
        or has_ste_vin_bnd
        or amchi_base.has_mobile_hydrogens(ich)
    )

    return is_bad


# # AMChI layer functions
# # # Formula layer
def formula_string(gra):
    """AMChI formula layer from graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    """
    fml = formula(gra)
    fml_str = form.string(fml)
    return fml_str


# # # Main layers
def connection_layer(gra, default_old: bool = True):
    """AMChI connection (c) layer from graph.

    To do: Remove the default_old flag and use only this implementation

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param default_old: Default to the old implementation, for non-polycycles?
    :returns: the connection layer, without prefix
    :rtype: str
    """
    # For now, use the old connection layer code for all but polycycles
    if default_old and len(rings_atom_keys(gra)) <= 1:
        return connection_layer_deprecated(gra)

    # Switch to one-indexing
    gra = relabel(gra, {k: k + 1 for k in atom_keys(gra)})

    # If there are terminal atoms, start from the one with the lowest canonical
    # number
    term_keys = terminal_atom_keys(gra, backbone=True)
    start_key = min(term_keys) if term_keys else 1

    # Identify children and parents for depth-first traversal
    dfs = dfs_(gra, start_key)
    child_dct = dfs_children(dfs)
    parent_dct = dfs_parents(dfs)

    # Identify ring bonds that won'tibe visited in the depth-first traversal
    # These will be connected back from the later DFS atom to the earlier one
    miss_child_dct = dfs_missing_children(gra, dfs)

    # Perform depth-first traversal, building the connection layer
    branch_depth = 0
    branch_start_keys = set()
    keys = [start_key]
    conn_lyr = ""
    while keys:
        key = keys.pop()

        # 1. Open parentheses for a branch, or replace ')' with ',' if adding additional
        # branches from the same parent
        if key in branch_start_keys:
            branch_depth += 1
            branch_start_keys.remove(key)
            if conn_lyr.endswith(")"):
                conn_lyr = re.sub("[)]$", ",", conn_lyr)
            else:
                conn_lyr += "("
        # 2. Encode the bond to the parent atom
        elif key in parent_dct and not conn_lyr.endswith(")"):
            conn_lyr += "-"

        # 3. Encode the atom
        conn_lyr += f"{key}"

        # 4. Encode missing ring closure bonds
        if key in miss_child_dct:
            rng_keys = miss_child_dct.pop(key)
            term_key = None if key in child_dct else rng_keys.pop()
            if rng_keys:
                conn_lyr += "(" + ",".join(map(str, rng_keys)) + ")"
            if term_key is not None:
                conn_lyr += f"{term_key}" if conn_lyr.endswith(")") else f"-{term_key}"

        # 5. Set continuation atoms, or terminate the string or branch
        if key in child_dct:
            child_keys = child_dct[key]
            branch_start_keys.update(child_keys[1:])
            keys.extend(child_keys)
        elif branch_depth:
            branch_depth -= 1
            conn_lyr += ")"

    conn_lyr = conn_lyr if conn_lyr != "1" else ""
    return conn_lyr


def connection_layer_deprecated(gra):
    """AMChI connection (c) layer from graph

    Deprecated: Needs to be dropped, but will require database changes

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: the connection layer, without prefix
    :rtype: str
    """
    # Get a one-indexed neighbor keys dictionary.
    nkeys_dct = {
        k + 1: [n + 1 for n in ns] for k, ns in atoms_neighbor_atom_keys(gra).items()
    }

    # Get a one-indexed list of ring keys
    rng_keys_lst = [[k + 1 for k in ks] for ks in rings_atom_keys(gra)]

    def _recurse_connection_layer(conn_lyr, conn_lst, key, just_seen=None):
        nkeys = nkeys_dct.pop(key) if key in nkeys_dct else []

        # Remove keys just seen from the list of neighbors, to avoid doubling
        # back.
        if just_seen in nkeys:
            nkeys.remove(just_seen)

        # Start the connection layer (string) and list. We could just work with
        # the layer directly, but the list is necessary for sorting.
        conn_lyr = f"{key}"
        conn_lst = [key]

        # Now, extend the layer/list along the neighboring atoms.
        if nkeys:
            # Build sub-layers/lists by recursively calling this function.
            sub_lyrs = []
            sub_lsts = []
            while nkeys:
                nkey = nkeys.pop(0)
                sub_lyr, sub_lst = _recurse_connection_layer(
                    "", [], nkey, just_seen=key
                )

                sub_lyrs.append(sub_lyr)
                sub_lsts.append(sub_lst)

                # If this is a ring, remove the neighbor on the other side of
                # `key` to prevent repetition as we go around the ring.
                # I can't think of a cleaner way to do this right now.
                flat_sub_lst = list(util.flatten(sub_lst))
                if key in flat_sub_lst:
                    for rng_keys in rng_keys_lst:
                        if key in rng_keys and nkey in rng_keys:
                            rkeys = util.ring.cycle_item_to_front(
                                rng_keys, key, end_item=nkey
                            )
                            if rkeys[1] in nkeys:
                                nkeys.remove(rkeys[1])

            # Now, join the sub-layers and lists together.
            # If there is only one neighbor, we join it as
            #   k-n-...
            if len(sub_lsts) == 1:
                # Extend the layer string
                conn_lyr += f"-{sub_lyrs[0]}"

                # Extend the list
                sub_lst = sub_lsts[0]
                conn_lst.extend(sub_lst)
            # If there are multiple neighbors, we join it as
            #   k(n1-...,n2-...)n3-...
            else:
                # Sort the list of branches by length and index values.
                srt_idxs = sorted(
                    range(len(sub_lsts)), key=lambda i: (len(sub_lsts[i]), sub_lsts[i])
                )

                # Apply the sort to both layers and lists.
                sub_lyrs = list(map(sub_lyrs.__getitem__, srt_idxs))
                sub_lsts = list(map(sub_lsts.__getitem__, srt_idxs))

                # Extend the layer string.
                conn_lyr += f"({','.join(sub_lyrs[:-1])}){sub_lyrs[-1]}"

                # Append the lists of neighboring branches.
                conn_lst.append(sub_lsts)

        return conn_lyr, conn_lst

    # If there are terminal atoms, start from the one with the lowest canonical
    # number
    term_keys = terminal_atom_keys(gra, backbone=True)
    start_key = min(term_keys) + 1 if term_keys else 1
    conn_lyr, _ = _recurse_connection_layer("", [], start_key)
    conn_lyr = conn_lyr if conn_lyr != "1" else ""

    return conn_lyr


def hydrogen_layer(gra):
    """AMChI hydrogen (h) layer from graph

    :param gra: implicit molecular graph
    :type gra: automol graph data structure
    :returns: the hydrogen layer, without prefix
    :rtype: str
    """
    # Determine hydrogen counts
    nhyd_dct = atom_implicit_hydrogens(gra)
    all_keys = sorted(atom_keys(gra), key=nhyd_dct.__getitem__)
    grps = [
        (nh, sorted(k + 1 for k in ks))
        for nh, ks in itertools.groupby(all_keys, key=nhyd_dct.__getitem__)
        if nh > 0
    ]

    # Build the hydrogen layer string
    slyrs = []
    for nhyd, keys in grps:
        parts = util.equivalence_partition(keys, lambda x, y: y in (x - 1, x + 1))
        parts = sorted(map(sorted, parts))
        strs = [
            "{:d}-{:d}".format(min(p), max(p)) if len(p) > 1 else "{:d}".format(p[0])
            for p in parts
        ]
        if nhyd == 1:
            slyrs.append(",".join(strs) + "H")
        else:
            slyrs.append(",".join(strs) + f"H{nhyd}")

    nhyd_lyr = ",".join(slyrs)
    return nhyd_lyr


# # # Stereo layers
def bond_stereo_layer(gra):
    """AMChI bond stereo (b) layer from graph

    cis     = '-' = False
    trans   = '+' = True

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: the bond stereo layer, without prefix
    :rtype: str
    """
    # Generate the bond stereo layer string
    bnd_keys = sorted(sorted(k, reverse=True) for k in bond_stereo_keys(gra))
    bnd_par_dct = bond_stereo_parities(gra)

    bnd_pars = list(map(bnd_par_dct.__getitem__, map(frozenset, bnd_keys)))
    bnd_sgns = ["+" if p else "-" for p in bnd_pars]
    bnd_strs = [f"{i+1}-{j+1}{s}" for (i, j), s in zip(bnd_keys, bnd_sgns)]

    bnd_ste_lyr = ",".join(bnd_strs)
    return bnd_ste_lyr


def atom_stereo_layer(gra):
    """AMChI atom stereo (t, m) layer from graph

    S = counterclockwise = '@'  = '-' = False
    R = clockwise        = '@@' = '+' = True

    :param gra: molecular graph
    :type gra: automol graph data structure
    :rtype: (str, str)
    """
    # Generate the atom stereo layer strings
    atm_keys = sorted(atom_stereo_keys(gra))
    atm_par_dct = atom_stereo_parities(gra)

    atm_pars = list(map(atm_par_dct.__getitem__, atm_keys))
    atm_sgns = ["+" if p else "-" for p in atm_pars]
    atm_strs = [f"{i+1}{s}" for i, s in zip(atm_keys, atm_sgns)]

    atm_ste_lyr = ",".join(atm_strs)
    return atm_ste_lyr


# # # TS layers
def bond_breaking_layer(gra):
    """AMChI bond breaking (k) layer from graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: the bond breaking layer, without prefix
    :rtype: str
    """
    # Generate the bond stereo layer string
    brk_bnd_keys = sorted(sorted(k, reverse=True) for k in ts_breaking_bond_keys(gra))

    brk_bnd_strs = [f"{i+1}-{j+1}" for (i, j) in brk_bnd_keys]

    brk_bnd_lyr = ",".join(brk_bnd_strs)
    return brk_bnd_lyr


def bond_forming_layer(gra):
    """AMChI bond forming (f) layer from graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: the bond forming layer, without prefix
    :rtype: str
    """
    # Generate the bond stereo layer string
    frm_bnd_keys = sorted(sorted(k, reverse=True) for k in ts_forming_bond_keys(gra))

    frm_bnd_strs = [f"{i+1}-{j+1}" for (i, j) in frm_bnd_keys]

    frm_bnd_lyr = ",".join(frm_bnd_strs)
    return frm_bnd_lyr
