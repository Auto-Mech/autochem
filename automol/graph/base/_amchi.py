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
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._core import terminal_heavy_atom_keys
from automol.graph.base._algo import is_connected
from automol.graph.base._canon import canonical_enantiomer
import automol.amchi.base


# AMChI functions
def amchi(gra, stereo=True, can=True, is_reflected=None):
    """ AMChI string from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param stereo: Include stereo in the AMChI string, if present?
        :type stereo: bool
        :param can: Canonicalize the graph? Set to True by default, causing the
            graph to be canonicalized. If setting to False to avoid
            re-canonicalization, the `is_reflected` flag must be set for a
            canonical result.
        :type can: bool
        :param is_reflected: If using pre-canonicalized graph, is it a
            reflected enantiomer? If True, yes; if False, it's an enantiomer
            that isn't reflected; if None, it's not an enantiomer.
        :type is_reflected: bool or NoneType
        :returns: the AMChI string
        :rtype: str
    """
    assert is_connected(gra), (
        "Cannot form connection layer for disconnected graph.")

    if not stereo:
        gra = without_stereo_parities(gra)

    # Convert to implicit graph
    gra = implicit(gra)

    # Canonicalize and determine canonical enantiomer
    if can:
        gra, is_reflected = canonical_enantiomer(gra)

    fml_str = _formula_string(gra)
    main_lyr_dct = _main_layers(gra)
    ste_lyr_dct = _stereo_layers(gra, is_reflected=is_reflected)

    chi = automol.amchi.base.from_data(fml_str=fml_str,
                                       main_lyr_dct=main_lyr_dct,
                                       ste_lyr_dct=ste_lyr_dct)
    return chi


# # AMChI layer functions
# # # Formula layer
def _formula_string(gra):
    """ AMChI formula layer from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
    """
    fml = formula(gra)
    fml_str = automol.formula.string(fml)
    return fml_str


# # # Main layers
def _main_layers(gra):
    """ Determine the main layers, describing the connectivity of the molecule.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: the 'c' and 'h' layers, as a dictionary
        :rtype: str
    """
    conn_lyr = _connection_layer(gra)
    nhyd_lyr = _hydrogen_layer(gra)
    lyr_dct = {'c': conn_lyr, 'h': nhyd_lyr}
    return lyr_dct


def _connection_layer(gra):
    """ AMChI connection (c) layer from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: the connection layer, without prefix
        :rtype: str
    """
    conn_lyr, _ = _connection_layer_and_list(gra)
    return conn_lyr


def _hydrogen_layer(gra):
    """ AMChI hydrogen (h) layer from graph

        :param gra: implicit molecular graph
        :type gra: automol graph data structure
        :returns: the hydrogen layer, without prefix
        :rtype: str
    """
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


def _connection_layer_and_list(gra):
    """ AMChI connection layer and list from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: the connection layer, without prefix, and connection list
        :rtype: str, list
    """
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

                # Append the lists of neighboring branches.
                conn_lst.append(sub_lsts)

        return conn_lyr, conn_lst

    # If there are terminal atoms, start from the one with the lowest canonical
    # number
    term_keys = terminal_heavy_atom_keys(gra)
    start_key = min(term_keys) + 1 if term_keys else 1
    conn_lyr, conn_lst = _recurse_connection_layer('', [], start_key)
    return conn_lyr, conn_lst


# # # Stereo layers
def _stereo_layers(gra, is_reflected=None):
    """ Determine the stereo layers, describing bond and atom stereochemistry.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param is_reflected: Is this a reflected enantiomer? If True, yes; if
            False, it's an enantiomer that isn't reflected; if None, it's not
            an enantiomer.
        :type is_reflected: bool or NoneType
        :returns: the 'b', 't', 'm', and 's' layers, as a dictionary
        :rtype: str
    """
    b_lyr = _bond_stereo_layer(gra)
    t_lyr = _atom_stereo_layer(gra)

    lyr_dct = {}
    if b_lyr:
        lyr_dct['b'] = b_lyr
    if t_lyr:
        lyr_dct['t'] = t_lyr
    if is_reflected is not None:
        assert t_lyr, (
            "If this is an enantiomer, there must be an atom stereo layer.")
        lyr_dct['m'] = '1' if is_reflected else '0'
        lyr_dct['s'] = '1'

    return lyr_dct


def _bond_stereo_layer(gra):
    """ AMChI bond stereo (b) layer from graph

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
    bnd_sgns = ['+' if p else '-' for p in bnd_pars]
    bnd_strs = [f'{i+1}-{j+1}{s}' for (i, j), s in zip(bnd_keys, bnd_sgns)]

    bnd_ste_lyr = ','.join(bnd_strs)
    return bnd_ste_lyr


def _atom_stereo_layer(gra):
    """ AMChI atom stereo (t, m) layer from graph

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
    atm_sgns = ['+' if p else '-' for p in atm_pars]
    atm_strs = [f'{i+1}{s}' for i, s in zip(atm_keys, atm_sgns)]

    atm_ste_lyr = ','.join(atm_strs)
    return atm_ste_lyr


if __name__ == '__main__':
    # bond stereo
    GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
            3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
           {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
            frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
            frozenset({2, 5}): (1, False)})
    CHI = amchi(GRA)
    print(CHI)
    ICH = automol.graph.inchi(GRA, stereo=True)
    print(ICH)

    # atom stereo
    GRA = ({0: ('C', 1, None), 1: ('C', 1, True), 2: ('C', 1, True),
            3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
            6: ('F', 0, None), 7: ('F', 0, None)},
           {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
            frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
            frozenset({2, 7}): (1, None)})
    CHI = amchi(GRA)
    print(CHI)
    ICH = automol.graph.inchi(GRA, stereo=True)
    print(ICH)
