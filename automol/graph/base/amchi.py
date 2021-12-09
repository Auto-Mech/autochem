""" AutoMechanic Chemical Identifier (AMChI) generating functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import implicit
from automol.graph.base._core import relabel
from automol.graph.base._canon import canonical_keys
from automol.graph.base._algo import is_connected


def connection_layer(gra, can_key_dct=None):
    """ AMChI connection layer from graph
    """
    assert is_connected(gra), (
        "Cannot form connection layer for disconnected graph.")

    def _join_with_dashes(nums):
        print(nums)
        return '-'.join(map('{:d}'.format, nums))

    # Don't recalculate canonical keys unless we have to
    can_key_dct = canonical_keys(gra) if can_key_dct is None else can_key_dct

    # Convert to a canonical graph
    gra = implicit(gra)
    gra = relabel(gra, can_key_dct)

    # Get a one-indexed neighbor keys dictionary
    nkeys_dct = {k+1: sorted(n+1 for n in ns) for k, ns in
                 atoms_neighbor_atom_keys(gra).items()}
    print(nkeys_dct)

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
                sub_lyrs = list(map(_join_with_dashes, sub_lsts))
                conn_lyr += f"({','.join(sub_lyrs[:-1])}){sub_lyrs[-1]}"

        return conn_lst, conn_lyr

    _, conn_lyr = _recurse_connection_list([], '', 1)
    return conn_lyr


if __name__ == '__main__':
    import automol

    ICH = automol.smiles.inchi('c1ccccc1C(CCC)(CCC)CCC')
    GRA = automol.inchi.graph(ICH)
    print(ICH)
    CONN_LYR = connection_layer(GRA)
    print(CONN_LYR)
