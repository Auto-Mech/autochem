""" AutoMechanic Chemical Identifier (AMChI) generating functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
from automol import util
import automol.formula
from automol.graph.base._core import atom_keys
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import formula
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import implicit
from automol.graph.base._core import explicit
from automol.graph.base._core import relabel
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._canon import canonical_keys
from automol.graph.base._algo import is_connected
from automol.graph.base._stereo import to_index_based_stereo
import automol.amchi.base


# # main functions
def amchi(gra, stereo=True, can_key_dct=None):
    """ AMChI string from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can_key_dct: optionally, pass in known canonical keys to avoid
            recalculating them; if None, they will be calculated
        :type can_key_dct: dict[int: int]
        :returns: the AMChI string
        :rtype: str
    """
    if stereo:
        raise NotImplementedError

    # Don't recalculate canonical keys unless we have to
    can_key_dct = canonical_keys(gra) if can_key_dct is None else can_key_dct

    fml_str = amchi_formula_string(gra)
    main_lyr_dct = amchi_main_layers(gra, can_key_dct=can_key_dct)

    ach = automol.amchi.base.from_data(fml_str=fml_str,
                                       main_lyr_dct=main_lyr_dct)

    return ach


# # AMChI layer functions
def amchi_formula_string(gra):
    """ AMChI formula layer from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
    """
    fml = formula(gra)
    fml_str = automol.formula.string(fml)
    return fml_str


def amchi_main_layers(gra, can_key_dct=None):
    """ Determine the main layers, describing the connectivity of the molecule.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can_key_dct: optionally, pass in known canonical keys to avoid
            recalculating them; if None, they will be calculated
        :type can_key_dct: dict[int: int]
        :returns: the 'c' and 'h' layers, as a dictionary
        :rtype: str
    """
    # Don't recalculate canonical keys unless we have to
    can_key_dct = canonical_keys(gra) if can_key_dct is None else can_key_dct

    conn_lyr = amchi_connection_layer(gra, can_key_dct=can_key_dct)
    nhyd_lyr = amchi_hydrogen_layer(gra, can_key_dct=can_key_dct)
    main_lyr_dct = {'c': conn_lyr, 'h': nhyd_lyr}
    return main_lyr_dct


def amchi_connection_layer(gra, can_key_dct=None):
    """ AMChI connection (c) layer from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can_key_dct: optionally, pass in known canonical keys to avoid
            recalculating them; if None, they will be calculated
        :type can_key_dct: dict[int: int]
        :returns: the connection layer, without prefix
        :rtype: str
    """
    conn_lyr, _ = _connection_layer_and_list(gra, can_key_dct=can_key_dct)
    return conn_lyr


def amchi_hydrogen_layer(gra, can_key_dct=None):
    """ AMChI hydrogen (h) layer from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can_key_dct: optionally, pass in known canonical keys to avoid
            recalculating them; if None, they will be calculated
        :type can_key_dct: dict[int: int]
        :returns: the hydrogen layer, without prefix
        :rtype: str
    """
    # Don't recalculate canonical keys unless we have to
    can_key_dct = canonical_keys(gra) if can_key_dct is None else can_key_dct

    # Convert to a canonical graph
    gra = implicit(gra)
    gra = relabel(gra, can_key_dct)

    # Determine hydrogen counts
    nhyd_dct = atom_implicit_hydrogen_valences(gra)
    all_keys = sorted(atom_keys(gra), key=nhyd_dct.__getitem__)
    grps = [(nh, sorted(k+1 for k in ks)) for nh, ks in
            itertools.groupby(all_keys, key=nhyd_dct.__getitem__) if nh > 0]

    # Build the hydrogen layer string
    slyrs = []
    for nhyd, keys in grps:
        parts = util.equivalence_partition(keys, lambda x, y: y in (x-1, x+1))
        parts = sorted(map(sorted, parts))
        strs = ['{:d}-{:d}'.format(min(p), max(p)) if len(p) > 1
                else '{:d}'.format(p[0]) for p in parts]
        if nhyd == 1:
            slyrs.append(','.join(strs) + 'H')
        else:
            slyrs.append(','.join(strs) + f'H{nhyd}')

    nhyd_lyr = ','.join(slyrs)
    return nhyd_lyr


def amchi_bond_stereo_layer(gra, can_key_dct=None, can_stereo=False):
    """ AMChI bond stereo (b) layer from graph

        trans   = '+' = False
        cis     = '-' = True

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can_key_dct: optionally, pass in known canonical keys to avoid
            recalculating them; if None, they will be calculated
        :type can_key_dct: dict[int: int]
        :param can_stereo: is this graph using canonical stereo? if not, we
            must convert to canonical stereo
        :type can_stereo: bool
        :returns: the bond stereo layer, without prefix
        :rtype: str
    """
    # Don't recalculate canonical keys unless we have to
    can_key_dct = canonical_keys(gra) if can_key_dct is None else can_key_dct

    # Convert to a canonical graph
    gra = implicit(gra)
    gra = relabel(gra, can_key_dct)

    # If stereo isn't canonical, convert to canonical stereo before proceeding.
    if not can_stereo:
        gra = explicit(gra)
        gra = to_index_based_stereo(gra)
        gra = implicit(gra)

    # Generate the bond stereo layer string
    bnd_keys = sorted(sorted(k, reverse=True) for k in bond_stereo_keys(gra))
    bnd_par_dct = bond_stereo_parities(gra)

    bnd_pars = list(map(bnd_par_dct.__getitem__, map(frozenset, bnd_keys)))
    bnd_sgns = ['-' if p else '+' for p in bnd_pars]
    bnd_strs = [f'{i+1}-{j+1}{s}' for (i, j), s in zip(bnd_keys, bnd_sgns)]

    bnd_ste_lyr = ','.join(bnd_strs)
    return bnd_ste_lyr


def amchi_atom_stereo_layers(gra, can_key_dct=None, can_stereo=False):
    """ AMChI atom stereo (t, m) layer from graph

        R = clockwise        = '@@' = '+' = False
        S = counterclockwise = '@'  = '-' = True

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can_key_dct: optionally, pass in known canonical keys to avoid
            recalculating them; if None, they will be calculated
        :type can_key_dct: dict[int: int]
        :param can_stereo: is this graph using canonical stereo? if not, we
            must convert to canonical stereo
        :type can_stereo: bool
        :returns: the atom stereo layers, t and m, without prefix
        :rtype: (str, str)
    """
    # Don't recalculate canonical keys unless we have to
    can_key_dct = canonical_keys(gra) if can_key_dct is None else can_key_dct

    # Convert to a canonical graph
    gra = implicit(gra)
    gra = relabel(gra, can_key_dct)

    # If stereo isn't canonical, convert to canonical stereo before proceeding.
    if not can_stereo:
        gra = explicit(gra)
        gra = to_index_based_stereo(gra)
        gra = implicit(gra)

    # Generate the atom stereo layer strings
    atm_keys = sorted(atom_stereo_keys(gra))
    atm_par_dct = atom_stereo_parities(gra)

    atm_pars = list(map(atm_par_dct.__getitem__, atm_keys))
    atm_sgns = ['-' if p else '+' for p in atm_pars]
    atm_strs = [f'{i+1}{s}' for i, s in zip(atm_keys, atm_sgns)]

    atm_ste_lyr = ','.join(atm_strs)
    return atm_ste_lyr


# # private helpers
def _connection_layer_and_list(gra, can_key_dct=None):
    """ AMChI connection layer and list from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can_key_dct: optionally, pass in known canonical keys to avoid
            recalculating them; if None, they will be calculated
        :type can_key_dct: dict[int: int]
        :returns: the connection layer, without prefix, and connection list
        :rtype: str, list
    """
    assert is_connected(gra), (
        "Cannot form connection layer for disconnected graph.")

    # Don't recalculate canonical keys unless we have to
    can_key_dct = canonical_keys(gra) if can_key_dct is None else can_key_dct

    # Convert to a canonical graph
    gra = implicit(gra)
    gra = relabel(gra, can_key_dct)

    # Get a one-indexed neighbor keys dictionary.
    nkeys_dct = {k+1: [n+1 for n in ns] for k, ns in
                 atoms_neighbor_atom_keys(gra).items()}

    def _recurse_connection_layer(conn_lyr, conn_lst, key, just_seen=None):
        nkeys = nkeys_dct.pop(key) if key in nkeys_dct else []

        # Remove keys just seen from the list of neighbors, to avoid doubling
        # back.
        if just_seen in nkeys:
            nkeys.remove(just_seen)

        # Start the connection layer (string) and list. We could just work with
        # the layer directly, but the list is necessary for sorting.
        conn_lyr = f'{key}'
        conn_lst = [key]

        # Now, extend the layer/list along the neighboring atoms.
        if nkeys:
            # Build sub-layers/lists by recursively calling this function.
            sub_lyrs = []
            sub_lsts = []
            while nkeys:
                nkey = nkeys.pop(0)
                sub_lyr, sub_lst = _recurse_connection_layer('', [], nkey,
                                                             just_seen=key)

                sub_lyrs.append(sub_lyr)
                sub_lsts.append(sub_lst)

                # If this is a ring, remove the neighbor on the other side of
                # `key` to prevent repetition as we go around the ring.
                if sub_lst[-1] == key:
                    nkeys.remove(sub_lst[-2])

            # Now, join the sub-layers and lists together.
            # If there is only one neighbor, we join it as
            #   k-n-...
            if len(sub_lsts) == 1:
                # Extend the layer string
                conn_lyr += f'-{sub_lyrs[0]}'

                # Extend the list
                sub_lst = sub_lsts[0]
                conn_lst.extend(sub_lst)
            # If there are multiple neighbors, we join it as
            #   k(n1-...,n2-...)n3-...
            else:
                # Sort the list of branches by length and index values.
                srt_idxs = sorted(
                    range(len(sub_lsts)),
                    key=lambda i: (len(sub_lsts[i]), sub_lsts[i]))

                # Apply the sort to both layers and lists.
                sub_lyrs = list(map(sub_lyrs.__getitem__, srt_idxs))
                sub_lsts = list(map(sub_lsts.__getitem__, srt_idxs))

                # Extend the layer string.
                conn_lyr += f"({','.join(sub_lyrs[:-1])}){sub_lyrs[-1]}"

                # Append the lists of neighboring brancches.
                conn_lst.append(sub_lsts)

        return conn_lyr, conn_lst

    conn_lyr, conn_lst = _recurse_connection_layer('', [], 1)
    return conn_lyr, conn_lst


if __name__ == '__main__':
    # # bond stereo
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
    #        {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
    #         frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
    #         frozenset({2, 5}): (1, False)})
    # CAN_KEY_DCT = dict(enumerate(range(6)))
    # print(amchi_bond_stereo_layer(GRA, can_key_dct=CAN_KEY_DCT))

    # atom stereo
    GRA = ({0: ('C', 1, None), 1: ('C', 1, False), 2: ('C', 1, False),
            3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
            6: ('F', 0, None), 7: ('F', 0, None)},
           {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
            frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
            frozenset({2, 7}): (1, None)})
    CAN_KEY_DCT = dict(enumerate(range(8)))
    amchi_atom_stereo_layers(GRA, can_key_dct=CAN_KEY_DCT)
