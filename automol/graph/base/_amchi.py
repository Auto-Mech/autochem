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
from automol.graph.base._core import relabel
from automol.graph.base._canon import canonical_keys
from automol.graph.base._algo import is_connected


def amchi_formula_layer(gra):
    """ AMChI formula layer from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
    """
    fml = formula(gra)
    fml_lyr = automol.formula.string(fml)
    return fml_lyr


def amchi_main_layers(gra):
    """ Determine the main layers, describing the connectivity of the molecule.

        :param gra: molecular graph
        :type gra: automol graph data structure
    """
    conn_lyr = amchi_connection_layer(gra)
    nhyd_lyr = amchi_hydrogen_layer(gra)
    main_lyr_dct = {'c': conn_lyr, 'h': nhyd_lyr}
    return main_lyr_dct


def amchi_connection_layer(gra, can_key_dct=None):
    """ AMChI connection layer from graph

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


def amchi_connection_list(gra, can_key_dct=None):
    """ AMChI connection list from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can_key_dct: optionally, pass in known canonical keys to avoid
            recalculating them; if None, they will be calculated
        :type can_key_dct: dict[int: int]
        :returns: the connection list
        :rtype: list
    """
    _, conn_lst = _connection_layer_and_list(gra, can_key_dct=can_key_dct)
    return conn_lst


def amchi_hydrogen_layer(gra, can_key_dct=None):
    """ AMChI hydrogen layer from graph

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
        slyrs.append(','.join(strs) + f'H{nhyd}')

    nhyd_lyr = ','.join(slyrs)
    return nhyd_lyr


# String generators
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

    def _join_with_dashes(nums):
        print('nums', nums)
        return '-'.join(map('{:d}'.format, nums))

    # Don't recalculate canonical keys unless we have to
    can_key_dct = canonical_keys(gra) if can_key_dct is None else can_key_dct

    # Convert to a canonical graph
    gra = implicit(gra)
    gra = relabel(gra, can_key_dct)

    # Get a one-indexed neighbor keys dictionary. Sort the neighbors so that
    # they come out in the right order when we form the connection layer.
    nkeys_dct = {k+1: sorted(n+1 for n in ns) for k, ns in
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
                # Extend the layer string.
                conn_lyr += f"({','.join(sub_lyrs[:-1])}){sub_lyrs[-1]}"

                # Append the lists of neighboring branches.  They have been
                # pre-sorted above where we formed the nkeys_dct, so for a
                # given length they will be sorted by index.
                sub_lsts = sorted(sub_lsts, key=len)
                conn_lst.append(sub_lsts)

        return conn_lyr, conn_lst

    conn_lyr, conn_lst = _recurse_connection_layer('', [], 1)
    return conn_lyr, conn_lst
