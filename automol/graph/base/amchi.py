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


def formula_layer(gra):
    """ AMChI formula layer from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
    """
    fml = formula(gra)
    fml_lyr = automol.formula.string(fml)
    return fml_lyr


def connection_layer(gra, can_key_dct=None):
    """ AMChI connection layer from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can_key_dct: optionally, pass in known canonical keys to avoid
            recalculating them; if None, they will be calculated
        :type can_key_dct: dict[int: int]
        :returns: the connection layer, without prefix
        :rtype: str
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

    # Get a one-indexed neighbor keys dictionary
    nkeys_dct = {k+1: sorted(n+1 for n in ns) for k, ns in
                 atoms_neighbor_atom_keys(gra).items()}

    def _recurse_connection_list(conn_lst, conn_lyr, key, just_seen=None):
        nkeys = nkeys_dct.pop(key) if key in nkeys_dct else []

        if just_seen in nkeys:
            nkeys.remove(just_seen)

        conn_lst = [key]
        conn_lyr = f'{key}'

        if nkeys:
            sub_lsts = []
            while nkeys:
                nkey = nkeys.pop(0)
                sub_lst, sub_lyr = _recurse_connection_list([], '', nkey,
                                                            just_seen=key)
                sub_lsts.append(sub_lst)

                # If this is a cycle, remove the neighbor on the other side of
                # `key` to prevent repetition
                if sub_lst[-1] == key:
                    nkeys.remove(sub_lst[-2])

            if len(sub_lsts) == 1:
                sub_lst = sub_lsts[0]
                conn_lst.extend(sub_lst)
                conn_lyr += f'-{sub_lyr}'
            else:
                sub_lsts = sorted(sub_lsts, key=len)
                conn_lst.append(sub_lsts)
                print(sub_lsts)
                # sub_lyrs = list(map(_join_with_dashes, sub_lsts))
                # conn_lyr += f"({','.join(sub_lyrs[:-1])}){sub_lyrs[-1]}"

        return conn_lst, conn_lyr

    conn_lst, conn_lyr = _recurse_connection_list([], '', 1)
    print(conn_lst)
    return conn_lyr


def connection_layer_old(gra, can_key_dct=None):
    """ AMChI connection layer from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can_key_dct: optionally, pass in known canonical keys to avoid
            recalculating them; if None, they will be calculated
        :type can_key_dct: dict[int: int]
        :returns: the connection layer, without prefix
        :rtype: str
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

    # Get a one-indexed neighbor keys dictionary
    nkeys_dct = {k+1: sorted(n+1 for n in ns) for k, ns in
                 atoms_neighbor_atom_keys(gra).items()}

    def _recurse_connection_list(conn_lst, conn_lyr, key, just_seen=None):
        nkeys = nkeys_dct.pop(key) if key in nkeys_dct else []

        if just_seen in nkeys:
            nkeys.remove(just_seen)

        conn_lst = [key]
        conn_lyr = f'{key}'

        if nkeys:
            sub_lsts = []
            while nkeys:
                nkey = nkeys.pop(0)
                sub_lst, sub_lyr = _recurse_connection_list([], '', nkey,
                                                            just_seen=key)
                sub_lsts.append(sub_lst)

                # If this is a cycle, remove the neighbor on the other side of
                # `key` to prevent repetition
                if sub_lst[-1] == key:
                    nkeys.remove(sub_lst[-2])

            if len(sub_lsts) == 1:
                sub_lst = sub_lsts[0]
                conn_lst.extend(sub_lst)
                conn_lyr += f'-{sub_lyr}'
            else:
                sub_lsts = sorted(sub_lsts, key=len)
                conn_lst.append(sub_lsts)
                print(sub_lsts)
                sub_lyrs = list(map(_join_with_dashes, sub_lsts))
                conn_lyr += f"({','.join(sub_lyrs[:-1])}){sub_lyrs[-1]}"

        return conn_lst, conn_lyr

    _, conn_lyr = _recurse_connection_list([], '', 1)
    return conn_lyr


def hydrogen_layer(gra, can_key_dct=None):
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


if __name__ == '__main__':
    import automol

    # ICH = automol.smiles.inchi('c1ccccc1C(CCCO)(CCCCl)CCC')
    ICH = automol.smiles.inchi('CC(C(C)C)C')
    GRA = automol.inchi.graph(ICH)
    print(ICH)
    print(formula_layer(GRA))
    print(connection_layer(GRA))
    print(hydrogen_layer(GRA))
