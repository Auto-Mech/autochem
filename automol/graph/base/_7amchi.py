""" AutoMechanic Chemical Identifier (AMChI) generating functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools

from automol import form, util
from automol.amchi import base as amchi_base
from automol.graph.base._0core import (
    atom_implicit_hydrogens,
    atom_keys,
    atom_stereo_keys,
    atom_stereo_parities,
    atoms_neighbor_atom_keys,
    bond_stereo_keys,
    bond_stereo_parities,
    formula,
    implicit,
    relabel,
    terminal_atom_keys,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    without_dummy_atoms,
    without_stereo,
)
from automol.graph.base._2algo import connected_components, rings_atom_keys
from automol.graph.base._3kekule import (
    kekules_bond_orders_collated,
    vinyl_radical_atom_keys,
)
from automol.graph.base._6canon import (
    break_priority_ties,
    canonical_enantiomer_with_keys,
    canonical_ts_direction,
)


# AMChI functions
def amchi(gra, stereo=True):
    """AMChI string from graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereo in the AMChI string, if present?
    :type stereo: bool
    :returns: the AMChI string
    :rtype: str
    """
    chi, _ = amchi_with_indices(gra, stereo=stereo)
    return chi


def amchi_with_indices(gra, stereo=True):
    """AMChI string and AMChI canonical indices from graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereo in the AMChI string, if present?
    :type stereo: bool
    :returns: the AMChI string and the AMChI canonical indices for each
        connected component (components in multi-component AMChI ordering)
    :rtype: (str, tuple[dct[int: int]])
    """
    gras = connected_components(gra)
    chis, chi_idx_dcts = zip(
        *(connected_amchi_with_indices(g, stereo=stereo) for g in gras)
    )
    srt_idxs = amchi_base.argsort(chis)
    chis = tuple(chis[i] for i in srt_idxs)
    chi_idx_dcts = tuple(chi_idx_dcts[i] for i in srt_idxs)
    chi = amchi_base.join(chis, sort=False)
    return chi, chi_idx_dcts


def connected_amchi_with_indices(gra, stereo=True, pri_dct=None, is_refl=None):
    """single-component AMChI string from a connected graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereo in the AMChI string, if present?
    :type stereo: bool
    :param pri_dct: Optionally, pass in canonical priorities to avoid
        recalculating. If this is an enantiomer, `gra` and `pri_dct` must
        reflect the *canonical enantiomer*, and the `is_refl` flag must be
        set.
    :type pri_dct: dict[int: int]
    :param is_refl: If using pre-canonicalized graph, is it a
        reflected enantiomer? If True, yes; if False, it's an enantiomer
        that isn't reflected; if None, it's not an enantiomer.
    :type is_refl: bool or NoneType
    :returns: the AMChI string
    :rtype: str
    """
    gra = without_dummy_atoms(gra)

    if not stereo:
        gra = without_stereo(gra)

    # Convert to implicit graph
    gra = implicit(gra)

    # 1. Identify the canonical direction for TS graphs (`None` for non-TS graphs)
    gra, is_rev = canonical_ts_direction(gra)

    # 2. Canonicalize and determine canonical enantiomer
    if pri_dct is None:
        gra, chi_idx_dct, is_refl = canonical_enantiomer_with_keys(gra)
    else:
        chi_idx_dct = break_priority_ties(gra, pri_dct)

    gra = relabel(gra, chi_idx_dct)

    # 3. Generate the appropriate layers
    fml_str = _formula_string(gra)
    main_lyr_dct = _main_layers(gra)
    ste_lyr_dct = _stereo_layers(gra, is_refl=is_refl)
    ts_lyr_dct = _ts_layers(gra, is_rev=is_rev)

    chi = amchi_base.from_data(
        fml_lyr=fml_str,
        main_lyr_dct=main_lyr_dct,
        ste_lyr_dct=ste_lyr_dct,
        ts_lyr_dct=ts_lyr_dct,
    )

    return chi, chi_idx_dct


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
def _formula_string(gra):
    """AMChI formula layer from graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    """
    fml = formula(gra)
    fml_str = form.string(fml)
    return fml_str


# # # Main layers
def _main_layers(gra):
    """Determine the main layers, describing the connectivity of the molecule.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: the 'c' and 'h' layers, as a dictionary
    :rtype: str
    """
    conn_lyr = _connection_layer(gra)
    nhyd_lyr = _hydrogen_layer(gra)
    lyr_dct = {"c": conn_lyr, "h": nhyd_lyr}
    return lyr_dct


def _connection_layer(gra):
    """AMChI connection (c) layer from graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: the connection layer, without prefix
    :rtype: str
    """
    conn_lyr, _ = _connection_layer_and_list(gra)
    return conn_lyr


def _hydrogen_layer(gra):
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


def _connection_layer_and_list(gra):
    """AMChI connection layer and list from graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: the connection layer, without prefix, and connection list
    :rtype: str, list
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
    conn_lyr, conn_lst = _recurse_connection_layer("", [], start_key)
    conn_lyr = conn_lyr if conn_lyr != "1" else ""

    return conn_lyr, conn_lst


# # # Stereo layers
def _stereo_layers(gra, is_refl=None):
    """Determine the stereo layers, describing bond and atom stereochemistry.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param is_refl: Is this a reflected enantiomer? If True, yes; if
        False, it's an enantiomer that isn't reflected; if None, it's not
        an enantiomer.
    :type is_refl: bool or NoneType
    :returns: the 'b', 't', 'm', and 's' layers, as a dictionary
    :rtype: str
    """
    b_lyr = _bond_stereo_layer(gra)
    t_lyr = _atom_stereo_layer(gra)

    lyr_dct = {}
    if b_lyr:
        lyr_dct["b"] = b_lyr
    if t_lyr:
        lyr_dct["t"] = t_lyr
    if is_refl is not None:
        assert t_lyr, "If this is an enantiomer, there must be an atom stereo layer."
        lyr_dct["m"] = "1" if is_refl else "0"
        lyr_dct["s"] = "1"

    return lyr_dct


def _bond_stereo_layer(gra):
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


def _atom_stereo_layer(gra):
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
def _ts_layers(gra, is_rev=None):
    """Determine the TS layers, describing breaking/forming bonds and TS direction

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param is_rev: Is this a reverse TS? If True, yes; if False, it's a TS that
        isn't reversed; if None, it's not a TS.
    :type is_rev: bool or NoneType
    :returns: the 'f', 'k', and 'r' layers, as a dictionary
    :rtype: dict
    """
    k_lyr = _bond_breaking_layer(gra)
    f_lyr = _bond_forming_layer(gra)

    lyr_dct = {}
    if k_lyr:
        lyr_dct["k"] = k_lyr
    if f_lyr:
        lyr_dct["f"] = f_lyr
    if is_rev is not None:
        assert (
            f_lyr or k_lyr
        ), "If this is a TS, there must be breaking or forming bond layers."
        lyr_dct["r"] = "1" if is_rev else "0"

    return lyr_dct


def _bond_breaking_layer(gra):
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


def _bond_forming_layer(gra):
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
