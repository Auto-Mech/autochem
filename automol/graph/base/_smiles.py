""" SMILES (Simplified Molecular Input Line Entry System) functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import numpy
from phydat import ptab
from automol import util
from automol.util import dict_
from automol.graph.base._core import atom_keys
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_orders
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import terminal_atom_keys
from automol.graph.base._core import string
from automol.graph.base._core import implicit
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._core import add_bonded_atom
from automol.graph.base._core import set_atom_implicit_hydrogen_valences
from automol.graph.base._algo import is_connected
from automol.graph.base._algo import rings_atom_keys
from automol.graph.base._algo import cycle_ring_atom_key_to_front
from automol.graph.base._resonance import dominant_resonance
from automol.graph.base._resonance import radical_atom_keys_from_resonance
from automol.graph.base._canon import canonical


ORGANIC_SUBSET = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']


def smiles(gra, stereo=True, local_stereo=False, res_stereo=False):
    """ SMILES string from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param stereo: Include stereo?
        :type stereo: bool
        :param local_stereo: Is the graph using local stereo assignments? That
            is, are they based on atom keys rather than canonical keys?
        :type local_stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :returns: the SMILES string
        :rtype: str
    """
    assert is_connected(gra), (
        "Cannot form connection layer for disconnected graph.")

    if not stereo:
        gra = without_stereo_parities(gra)

    # If not using local stereo assignments, canonicalize the graph first.
    # From this point on, the stereo parities can be assumed to correspond to
    # the neighboring atom keys.
    if not local_stereo:
        gra = canonical(gra)

    # Convert to implicit graph
    gra = implicit(gra)

    # Insert hydrogens necessary for bond stereo
    gra = _insert_stereo_hydrogens(gra)

    # Find a dominant resonance
    rgr = dominant_resonance(gra)

    # Determine atom symbols
    symb_dct = atom_symbols(rgr)

    # Determine atom implicit hydrogens
    nhyd_dct = atom_implicit_hydrogen_valences(rgr)

    # Determine bond orders for this resonance
    bnd_ord_dct = bond_orders(rgr)

    # Find radical sites for this resonance
    rad_atm_keys = radical_atom_keys_from_resonance(rgr)

    # Determine neighbors
    nkeys_dct = atoms_neighbor_atom_keys(rgr)

    # Find stereo parities
    atm_par_dct = dict_.filter_by_value(
        atom_stereo_parities(rgr), lambda x: x is not None)
    bnd_par_dct = dict_.filter_by_value(
        bond_stereo_parities(rgr), lambda x: x is not None)

    # Remove stereo parities if requested
    if not res_stereo:
        print('before')
        print(bnd_par_dct)
        bnd_par_dct = dict_.filter_by_key(
            bnd_par_dct, lambda x: bnd_ord_dct[x] == 2)
        print('after')
        print(bnd_par_dct)
    else:
        raise NotImplementedError("Not yet implemented!")

    def _atom_representation(key, just_seen=None, nkeys=(), closures=()):
        symb = ptab.to_symbol(symb_dct[key])
        nhyd = nhyd_dct[key]

        needs_brackets = key in rad_atm_keys or symb not in ORGANIC_SUBSET

        hyd_rep = f'H{nhyd}' if nhyd > 1 else ('H' if nhyd == 1 else '')
        par_rep = ''

        if key in atm_par_dct:
            needs_brackets = True

            skeys = [just_seen]
            if nhyd:
                assert nhyd == 1
                skeys.append(-numpy.inf)
            if closures:
                skeys.extend(closures)
            skeys.extend(nkeys)

            can_par = atm_par_dct[key]
            smi_par = can_par ^ util.is_odd_permutation(skeys, sorted(skeys))
            par_rep = '@@' if smi_par else '@'

        if needs_brackets:
            rep = f'[{symb}{par_rep}{hyd_rep}]'
        else:
            rep = f'{symb}'

        return rep

    # Get the pool of stereo bonds for the graph and set up a dictionary for
    # storing the ending representation.
    ste_bnd_key_pool = list(bnd_par_dct.keys())
    drep_dct = {}

    def _bond_representation(key, just_seen=None):
        key0 = just_seen
        key1 = key

        # First, handle the bond order
        if key0 is None or key1 is None:
            rep = ''
        else:
            bnd_ord = bnd_ord_dct[frozenset({key0, key1})]
            if bnd_ord == 1:
                rep = ''
            elif bnd_ord == 2:
                rep = '='
            elif bnd_ord == 3:
                rep = '#'
            else:
                raise ValueError("Bond orders greater than 3 not permitted.")

        drep = drep_dct[(key0, key1)] if (key0, key1) in drep_dct else ''

        bnd_key = next((b for b in ste_bnd_key_pool if key1 in b), None)
        if bnd_key is not None:
            # We've encountered a new stereo bond, so remove it from the pool
            ste_bnd_key_pool.remove(bnd_key)

            # Determine the atoms involved
            key2, = bnd_key - {key1}
            nkey1s = set(nkeys_dct[key1]) - {key2}
            nkey2s = set(nkeys_dct[key2]) - {key1}

            nmax1 = max(nkey1s)
            nmax2 = max(nkey2s)

            nkey1 = just_seen if just_seen in nkey1s else nmax1
            nkey2 = nmax2

            # Determine parity
            can_par = bnd_par_dct[bnd_key]
            smi_par = can_par if nkey1 == nmax1 else not can_par

            # Determine bond directions
            drep1 = drep if drep else '/'
            if just_seen in nkey1s:
                drep = drep1
                flip = not smi_par
            else:
                drep_dct[(key1, nkey1)] = drep1
                flip = smi_par

            drep2 = _flip_direction(drep1, flip=flip)

            drep_dct[(key2, nkey2)] = drep2

        rep += drep

        # Second, handle directionality (bond stereo)
        return rep

    # Get the pool of rings for the graph and set up a dictionary for storing
    # their tags. As the SMILES is built, each next ring that is encountered
    # will be given a tag, removed from the pool, and transferred to the tag
    # dictionary.
    rng_pool = list(rings_atom_keys(rgr))
    rng_tag_dct = {}

    def _ring_representation_with_nkeys_and_closures(key, nkeys=()):
        nkeys = nkeys.copy()

        # Check for new rings in the ring pool. If a new ring is found, create
        # a tag, add it to the tags dictionary, and drop it from the rings
        # pool.
        for new_rng in rng_pool:
            if key in new_rng:
                # Choose a neighbor key for SMILES ring closure
                clos_nkey = sorted(set(new_rng) & set(nkeys))[0]

                # Add it to the ring tag dictionary with the current key first
                # and the closure key last
                tag = max(rng_tag_dct.values(), default=0) + 1
                assert tag < 10, (
                    f"Ring tag exceeds 10 for this graph:\n{string(gra)}")
                rng = cycle_ring_atom_key_to_front(new_rng, key, clos_nkey)
                rng_tag_dct[rng] = tag

                # Remove it from the pool of unseen rings
                rng_pool.remove(new_rng)

        tags = []
        closures = []
        for rng, tag in rng_tag_dct.items():
            if key == rng[-1]:
                nkeys.remove(rng[0])
                closures.append(rng[0])
                # Handle the special case where the last ring bond has stereo
                if (rng[-1], rng[0]) in drep_dct:
                    drep = drep_dct[(rng[-1], rng[0])]
                    tags.append(f'{drep}{tag}')
                else:
                    tags.append(f'{tag}')
            if key == rng[0]:
                nkeys.remove(rng[-1])
                closures.append(rng[-1])
                tags.append(f'{tag}')

        rrep = ''.join(map(str, tags))
        return rrep, nkeys, closures

    # Determine neighboring keys
    nkeys_dct_pool = dict_.transform_values(
        atoms_neighbor_atom_keys(rgr), sorted)

    def _recurse_smiles(smi, lst, key, just_seen=None):
        nkeys = nkeys_dct_pool.pop(key) if key in nkeys_dct_pool else []

        # Remove keys just seen from the list of neighbors, to avoid doubling
        # back.
        if just_seen in nkeys:
            nkeys.remove(just_seen)

        # Start the SMILES string and connection list. The connection list is
        # used for sorting.
        rrep, nkeys, closures = _ring_representation_with_nkeys_and_closures(
            key, nkeys)
        arep = _atom_representation(key, just_seen, nkeys, closures=closures)
        brep = _bond_representation(key, just_seen)
        smi = f'{brep}{arep}{rrep}'
        lst = [key]

        # Now, extend the layer/list along the neighboring atoms.
        if nkeys:
            # Build sub-strings/lists by recursively calling this function.
            sub_smis = []
            sub_lsts = []
            while nkeys:
                nkey = nkeys.pop(0)
                sub_smi, sub_lst = _recurse_smiles('', [], nkey, just_seen=key)

                sub_smis.append(sub_smi)
                sub_lsts.append(sub_lst)

                # If this is a ring, remove the neighbor on the other side of
                # `key` to prevent repetition as we go around the ring.
                if sub_lst[-1] == key:
                    nkeys.remove(sub_lst[-2])

            # Now, join the sub-layers and lists together.
            # If there is only one neighbor, we joint it as
            #   {arep1}{brep2}{arep2}...
            if len(sub_lsts) == 1:
                sub_smi = sub_smis[0]
                sub_lst = sub_lsts[0]

                # Extend the SMILES string
                smi += f'{sub_smi}'

                # Extend the list
                lst.extend(sub_lst)
            # If there are multiple neighbors, we joint them as
            #   {arep1}({brep2}{arep2}...)({brep3}{arep3}...){brep4}{arep4}...
            else:
                assert len(sub_lsts) > 1

                # Extend the SMILES string
                smi += (''.join(map("({:s})".format, sub_smis[:-1]))
                        + sub_smis[-1])

                # Append the lists of neighboring branches.
                lst.append(sub_lsts)

        return smi, lst

    # If there are terminal atoms, start from the first one
    atm_keys = atom_keys(rgr)
    term_keys = terminal_atom_keys(gra, heavy=False)
    start_key = min(term_keys) if term_keys else min(atm_keys)

    smi, _ = _recurse_smiles('', [], start_key)

    return smi


# helpers
def _insert_stereo_hydrogens(gra):
    """ Insert hydrogens necessary for bond stereo into an implicit graph.
        Hydrogens are given negative keys for proper stereo sorting
    """
    bnd_keys = bond_stereo_keys(gra)
    nkeys_dct = atoms_neighbor_atom_keys(gra)
    nhyd_dct = atom_implicit_hydrogen_valences(gra)
    next_key = -max(atom_keys(gra)) - 1
    for bnd_key in bnd_keys:
        key1, key2 = bnd_key
        nkey1s = nkeys_dct[key1] - {key2}
        nkey2s = nkeys_dct[key2] - {key1}
        for key, nkeys in [(key1, nkey1s), (key2, nkey2s)]:
            if not nkeys:
                assert nhyd_dct[key] == 1
                gra = add_bonded_atom(gra, 'H', key, next_key)
                gra = set_atom_implicit_hydrogen_valences(gra, {key: 0})

                next_key = next_key - 1

    return gra


def _flip_direction(drep, flip=True):
    """ Flip the direction of a directional bond representation

        :param drep: the directional bond representation, '/' or '\\'
        :type drep: str
        :param flip: Flip the representation? If False, don't flip it.
        :type flip: bool
        :returns: the new representation, flipped if requested
        :rtype: str
    """
    if not flip:
        ret = drep
    else:
        ret = '\\' if drep == '/' else '/'
    return ret


if __name__ == '__main__':
    import automol

    # # Rings:
    # ICH = automol.smiles.inchi('C123C(O1)(CCO2)CNC3')
    # ICH = automol.smiles.inchi('C1=C(O[O])NC1=C[CH2]')
    # ICH = automol.smiles.inchi('C1CC1')

    # Test stereo atoms:
    # ICH = automol.smiles.inchi('N[C@](C)(F)C(=O)O')
    # ICH = automol.smiles.inchi('N[C@H](C)C(=O)O')
    # ICH = automol.smiles.inchi('C[C@H]1CCCCO1')
    # GRA = automol.inchi.graph(ICH)

    # Test stereo bonds:
    # ICH = automol.smiles.inchi(r'CC/C=C/C=C/[CH]F')
    # ICH = automol.smiles.inchi(r'CC/C=C/C=C/CF')
    # ICH = automol.smiles.inchi(r'C1CCCCCCCCCC/N=N/1')
    # ICH = automol.smiles.inchi(r'[H]/N=N/[H]')
    # ICH = automol.smiles.inchi(r'[H]/N=N/N=N\[H]')
    ICH = automol.smiles.inchi(r'F[CH2]/C=C/F')
    GEO = automol.inchi.geometry(ICH)
    GRA = automol.geom.graph(GEO)
    SMI = smiles(GRA)
    print(SMI)

    ICH_SMI = automol.inchi.smiles(ICH)
    print(ICH)
    print(ICH_SMI)

    SICH = automol.smiles.inchi(SMI)
    print(SICH)
    # print(automol.graph.string(GRA, one_indexed=True))
    print(automol.graph.rings_atom_keys(GRA))
    assert SICH == ICH
