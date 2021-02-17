""" functions for working with rotational bonds and groups
"""
import itertools
from automol.graph._graph_base import atom_keys
from automol.graph._graph_base import bond_keys
from automol.graph._graph_base import atom_symbols
from automol.graph._graph_base import atom_implicit_hydrogen_valences
from automol.graph._graph import explicit
from automol.graph._graph import implicit
from automol.graph._graph import without_dummy_atoms
from automol.graph._graph import subgraph
from automol.graph._graph import connected_components
from automol.graph._graph import terminal_heavy_atom_keys
from automol.graph._graph import branch_atom_keys
from automol.graph._graph import atoms_neighbor_atom_keys
from automol.graph._graph import dummy_atoms_neighbor_atom_key
from automol.graph._graph import atom_neighbor_atom_key
from automol.graph._ring import rings_bond_keys
from automol.graph._res import resonance_dominant_bond_orders


def rotational_bond_keys(gra, lin_keys=None, with_h_rotors=True):
    """ get all rotational bonds for a graph

    :param gra: the graph
    :param lin_keys: keys to linear atoms in the graph
    """
    gra = explicit(gra)
    sym_dct = atom_symbols(gra)
    ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    bnd_ord_dct = resonance_dominant_bond_orders(gra)
    rng_bnd_keys = list(itertools.chain(*rings_bond_keys(gra)))

    def _is_rotational_bond(bnd_key):
        ngb_keys_lst = [ngb_keys_dct[k] - bnd_key for k in bnd_key]

        is_single = max(bnd_ord_dct[bnd_key]) <= 1
        has_neighbors = all(ngb_keys_lst)
        not_in_ring = bnd_key not in rng_bnd_keys

        is_h_rotor = any(set(map(sym_dct.__getitem__, ks)) == {'H'}
                         for ks in ngb_keys_lst)

        return is_single and has_neighbors and not_in_ring and (
            not is_h_rotor or with_h_rotors)

    rot_bnd_keys = frozenset(filter(_is_rotational_bond, bond_keys(gra)))

    lin_keys_lst = linear_segments_atom_keys(gra, lin_keys=lin_keys)
    dum_keys = tuple(atom_keys(gra, sym='X'))
    for keys in lin_keys_lst:
        bnd_keys = sorted((k for k in rot_bnd_keys if k & set(keys)),
                          key=sorted)

        # Check whether there are neighboring atoms on either side of the
        # linear segment
        excl_keys = set(keys) | set(dum_keys)

        end_key1 = atom_neighbor_atom_key(
            gra, keys[0], excl_atm_keys=excl_keys)

        excl_keys |= {end_key1}
        end_key2 = atom_neighbor_atom_key(
            gra, keys[-1], excl_atm_keys=excl_keys)

        end_keys = {end_key1, end_key2}
        ngb_keys_lst = [ngb_keys_dct[k] - excl_keys for k in end_keys]
        has_neighbors = all(ngb_keys_lst)

        if not has_neighbors:
            rot_bnd_keys -= set(bnd_keys)
        else:
            rot_bnd_keys -= set(bnd_keys[:-1])

    return rot_bnd_keys


def rotational_groups(gra, key1, key2, dummy=False):
    """ get the rotational groups for a given rotational axis

    :param gra: the graph
    :param key1: the first atom key
    :param key2: the second atom key
    """

    if not dummy:
        gra = without_dummy_atoms(gra)

    bnd_key = frozenset({key1, key2})
    assert bnd_key in bond_keys(gra)
    grp1 = branch_atom_keys(gra, key2, bnd_key) - {key1}
    grp2 = branch_atom_keys(gra, key1, bnd_key) - {key2}
    grp1 = tuple(sorted(grp1))
    grp2 = tuple(sorted(grp2))
    return grp1, grp2


def rotational_symmetry_number(gra, key1, key2, lin_keys=None):
    """ get the rotational symmetry number along a given rotational axis

    :param gra: the graph
    :param key1: the first atom key
    :param key2: the second atom key
    """
    ngb_keys_dct = atoms_neighbor_atom_keys(without_dummy_atoms(gra))
    imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(implicit(gra))

    axis_keys = {key1, key2}
    # If the keys are part of a linear chain, use the ends of that for the
    # symmetry number calculation
    lin_keys_lst = linear_segments_atom_keys(gra, lin_keys=lin_keys)
    for keys in lin_keys_lst:
        if key1 in keys or key2 in keys:
            if len(keys) == 1:
                key1, key2 = sorted(ngb_keys_dct[keys[0]])
            else:
                key1, = ngb_keys_dct[keys[0]] - {keys[1]}
                key2, = ngb_keys_dct[keys[-1]] - {keys[-2]}
                axis_keys |= set(keys)
                break

    sym_num = 1
    for key in (key1, key2):
        if key in imp_hyd_vlc_dct:
            ngb_keys = ngb_keys_dct[key] - axis_keys
            if len(ngb_keys) == imp_hyd_vlc_dct[key] == 3:
                sym_num = 3
                break
    return sym_num


def linear_segments_atom_keys(gra, lin_keys=None):
    """ atom keys for linear segments in the graph
    """
    ngb_keys_dct = atoms_neighbor_atom_keys(without_dummy_atoms(gra))

    lin_keys = (dummy_atoms_neighbor_atom_key(gra).values()
                if lin_keys is None else lin_keys)

    lin_keys = [k for k in lin_keys if len(ngb_keys_dct[k]) <= 2]

    lin_segs = connected_components(subgraph(gra, lin_keys))

    lin_keys_lst = []
    for lin_seg in lin_segs:
        lin_seg_keys = atom_keys(lin_seg)
        if len(lin_seg_keys) == 1:
            key, = lin_seg_keys
            lin_keys_lst.append([key])
        else:
            end_key1, end_key2 = sorted(terminal_heavy_atom_keys(lin_seg))
            ngb_keys_dct = atoms_neighbor_atom_keys(lin_seg)

            key = None
            keys = [end_key1]
            while key != end_key2:
                key, = ngb_keys_dct[keys[-1]] - set(keys)
                keys.append(key)
            lin_keys_lst.append(keys)

    lin_keys_lst = tuple(map(tuple, lin_keys_lst))
    return lin_keys_lst


if __name__ == '__main__':
    import automol
    TSG = ({0: ('C', 0, None), 1: ('H', 0, None), 2: ('H', 0, None),
            3: ('H', 0, None), 4: ('H', 0, None), 5: ('X', 0, None),
            6: ('H', 0, None)},
           {frozenset({4, 5}): (0, None), frozenset({0, 3}): (1, None),
            frozenset({4, 6}): (0.1, None), frozenset({0, 1}): (1, None),
            frozenset({0, 2}): (1, None), frozenset({0, 4}): (0.9, None)})
    print(automol.graph.string(TSG, one_indexed=False))
    LIN_KEYS = [4]
    BND_KEYS = rotational_bond_keys(TSG, lin_keys=LIN_KEYS)
    print(BND_KEYS)
