""" graph-based z-matrix builder
"""
import automol.vmat
from automol.graph.base import atom_keys
from automol.graph.base import atom_symbols
from automol.graph.base import remove_bonds
from automol.graph.base import atoms_neighbor_atom_keys
from automol.graph.base import atoms_sorted_neighbor_atom_keys


def ring(gra, keys):
    """ generate a z-matrix for a ring

    All neighboring atoms along the ring will be included

    :param gra: the graph
    :param keys: ring keys, in the order they should appear in the z-matrix
    """
    # Break the bond between the first and last atoms to make this a chain
    gra = remove_bonds(gra, [(keys[0], keys[-1])])

    # Now, construct a v-matrix for the chain
    vma, zma_keys = chain(gra, keys, term_hydrogens=True)
    return vma, zma_keys


def chain(gra, keys, term_hydrogens=True):
    """ generate a v-matrix for a chain

    All neighboring atoms along the chain will be included

    :param gra: the graph
    :param keys: a list of keys for the chain
    :param term_hydrogens: whether or not to extend the chain to include
        terminal hydrogens, if present
    """
    if term_hydrogens:
        keys = _extend_chain_to_include_terminal_hydrogens(gra, keys)

    # 1. Start the chain
    vma, zma_keys = start_at(gra, keys[0])

    start_key, = set(keys) & set(_atoms_missing_neighbors(gra, zma_keys))
    keys = keys[keys.index(start_key):]

    # 2. Continue the chain
    if keys:
        vma, zma_keys = continue_chain(gra, keys, vma=vma, zma_keys=zma_keys)

    return vma, zma_keys


def continue_chain(gra, keys, vma, zma_keys, term_hydrogens=True):
    """ continue constructing a v-matrix along a chain

    All neighboring atoms along the chain will be included

    Exactly one atom in the chain must already be in the v-matrix

    :param gra: the graph for which the v-matrix will be constructed
    :param keys: the keys for atoms along the chain, which must be contiguous;
        the first atom must already appear in the v-matrix
    :param vma: a partial v-matrix from which to continue
    :param zma_keys: row keys for the partial v-matrix, identifying the atom
        specified by each row of `vma` in order
    :param extend: whether to extend the chain's start to include the three
        anchoring atoms
    :param term_hydrogens: whether to extend the chain's end to include
        terminal hydrogens
    """
    if term_hydrogens:
        keys = _extend_chain_to_include_terminal_hydrogens(gra, keys,
                                                           start=False)

    vma, zma_keys = complete_branch(gra, keys[0], vma, zma_keys,
                                    branch_keys=keys)

    return vma, zma_keys


def start_at(gra, key):
    """ start a v-matrix at a specific atom

    Returns the started vmatrix, along with keys to atoms whose neighbors are
    missing from it
    """
    symb_dct = atom_symbols(gra)
    ngb_keys_dct = atoms_sorted_neighbor_atom_keys(
        gra, symbs_first=('X', 'C',), symbs_last=('H',), ords_last=(0.1,))

    ngb_keys = ngb_keys_dct[key]
    if not ngb_keys:
        zma_keys = []
    elif len(ngb_keys) == 1:
        # Need special handling for atoms with only one neighbor
        if symb_dct[key] in ('H', 'X'):
            key2 = ngb_keys[0]
            zma_keys = (key2,) + ngb_keys_dct[key2]
        else:
            key2 = ngb_keys[0]
            ngb_keys = tuple(k for k in ngb_keys_dct[key2] if k != key)
            zma_keys = (key, key2) + ngb_keys
    else:
        zma_keys = (key,) + ngb_keys_dct[key]

    vma = ()
    for row, key_ in enumerate(zma_keys):
        idx1 = idx2 = idx3 = None
        if row > 0:
            key1 = next(k for k in ngb_keys_dct[key_] if k in zma_keys[:row])
            idx1 = zma_keys.index(key1)
        if row > 1:
            key2 = next(k for k in ngb_keys_dct[key1] if k in zma_keys[:row]
                        and k != key_)
            idx2 = zma_keys.index(key2)
        if row > 2:
            key3 = next(k for k in zma_keys[:row]
                        if k not in (key_, key1, key2))
            idx3 = zma_keys.index(key3)

        sym = symb_dct[key_]
        key_row = [idx1, idx2, idx3]
        vma = automol.vmat.add_atom(vma, sym, key_row)

    return vma, zma_keys


def complete_branch(gra, key, vma, zma_keys, branch_keys=None):
    """ continue constructing a v-matrix along a chain

    All neighboring atoms along the chain will be included

    Exactly one atom in the chain must already be in the v-matrix

    :param gra: the graph for which the v-matrix will be constructed
    :param vma: a partial v-matrix from which to continue
    :param zma_keys: row keys for the partial v-matrix, identifying the atom
        specified by each row of `vma` in order
    :param branch_keys: optionally, restrict the v-matrix to these keys and
        their neighbors; if `None`, the entire branch will be included
    """
    branch_keys = atom_keys(gra) if branch_keys is None else branch_keys
    keys = _extend_chain_to_include_anchoring_atoms(gra, [key], zma_keys)

    zma_keys = list(zma_keys)
    symb_dct = atom_symbols(gra)
    ngb_keys_dct = atoms_sorted_neighbor_atom_keys(
        gra, symbs_first=('X', 'C',), symbs_last=('H',), ords_last=(0.1,))

    def _continue(key1, key2, key3, vma, zma_keys):
        k3ns = list(ngb_keys_dct[key3])
        for k3n in set(k3ns) & set(zma_keys):
            k3ns.remove(k3n)

        if k3ns:
            key4 = k3ns.pop(0)

            # Add the leading atom to the v-matrix
            symb = symb_dct[key4]
            key_row = list(map(zma_keys.index, (key3, key2, key1)))
            vma = automol.vmat.add_atom(vma, symb, key_row)
            assert key4 not in zma_keys, f"Atom {key4:d} already in v-matrix"
            zma_keys.append(key4)

            # Add the neighbors of atom 3 (if any) to the v-matrix, decoupled
            # from atom 1 for properly decopuled torsions
            for k3n in k3ns:
                sym = symb_dct[k3n]

                if symb_dct[key4] == 'X':
                    key_row = list(map(zma_keys.index, (key3, key4, key2)))
                else:
                    key_row = list(map(zma_keys.index, (key3, key2, key4)))

                vma = automol.vmat.add_atom(vma, sym, key_row)
                assert k3n not in zma_keys, f"Atom {k3n:d} already in v-matrix"
                zma_keys.append(k3n)

            # Recursion
            if key4 in branch_keys:
                vma, zma_keys = _continue(key2, key3, key4, vma, zma_keys)

            if symb_dct[key4] == 'X':
                key2 = key4

            for k3n in k3ns:
                if k3n in branch_keys:
                    vma, zma_keys = _continue(key2, key3, k3n, vma, zma_keys)

        return vma, zma_keys

    key1, key2, key3 = keys[0], keys[1], keys[2]
    vma, zma_keys = _continue(key1, key2, key3, vma, zma_keys)

    return vma, zma_keys


# helpers
def _extend_chain_to_include_anchoring_atoms(gra, keys, zma_keys):
    """ extend chain to include three atoms already specified in v-matrix

    :param gra: the graph
    :param keys: keys in the chain; the first atom should already be specified
    :param zma_keys: keys currently in the v-matrix
    """
    ngb_keys_dct = atoms_sorted_neighbor_atom_keys(
        gra, symbs_first=('X', 'C',), symbs_last=('H',), ords_last=(0.1,))

    key3 = keys[0]
    assert key3 in zma_keys
    key2 = next(k for k in ngb_keys_dct[key3] if k in zma_keys)
    key1 = next(k for k in ngb_keys_dct[key2] if k in zma_keys and k != key3)
    keys = (key1, key2,) + tuple(keys)

    return keys


def _extend_chain_to_include_terminal_hydrogens(gra, keys,
                                                start=True, end=True):
    """ extend each end of a chain to include terminal hydrogens, if any
    """
    symb_dct = atom_symbols(gra)
    atm_ngb_dct = atoms_neighbor_atom_keys(gra)

    sta_ngbs = atm_ngb_dct[keys[0]] - {keys[1]}
    end_ngbs = atm_ngb_dct[keys[-1]] - {keys[-2]}

    sta_ngb = min((k for k in sta_ngbs if symb_dct[k] == 'H'), default=None)
    end_ngb = min((k for k in end_ngbs if symb_dct[k] == 'H'), default=None)

    keys = tuple(keys)

    if start and sta_ngb is not None:
        keys = (sta_ngb,) + keys

    if end and end_ngb is not None:
        keys = keys + (end_ngb,)

    return keys


def _atoms_missing_neighbors(gra, zma_keys):
    """ get atoms from the list currently in the v-matrix with neighbors that
    are not in the v-matrix
    """
    ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    keys = []
    for key in zma_keys:
        if any(k not in zma_keys for k in ngb_keys_dct[key]):
            keys.append(key)
    keys = tuple(keys)
    return keys
