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
# from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import terminal_heavy_atom_keys
from automol.graph.base._core import string
from automol.graph.base._core import implicit
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._algo import is_connected
from automol.graph.base._algo import rings_atom_keys
from automol.graph.base._algo import cycle_ring_atom_key_to_front
from automol.graph.base._resonance import dominant_resonance
from automol.graph.base._resonance import radical_atom_keys_from_resonance
from automol.graph.base._canon import canonical


ORGANIC_SUBSET = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']


def smiles(gra, stereo=True, local_stereo=False):
    """ SMILES string from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param stereo: Include stereo?
        :type stereo: bool
        :param local_stereo: Is the graph using local stereo assignments? That
            is, are they based on atom keys rather than canonical keys?
        :type local_stereo: bool
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

    # Find stereo parities
    atm_par_dct = dict_.filter_by_value(
        atom_stereo_parities(rgr), lambda x: x is not None)
    # bnd_par_dct = dict_.filter_by_value(
    #     bond_stereo_parities(rgr), lambda x: x is not None)
    print(f'atm_par_dct {atm_par_dct}')

    def _atom_representation(key, just_seen=None, nkeys=(), closures=()):
        symb = ptab.to_symbol(symb_dct[key])
        nhyd = nhyd_dct[key]
        if key in atm_par_dct:
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
            rep = f'[{symb}{par_rep}H]' if nhyd else f'[{symb}{par_rep}]'
        elif key in rad_atm_keys or symb not in ORGANIC_SUBSET:
            rep = f'[{symb}H{nhyd}]' if nhyd else f'[{symb}]'
        else:
            rep = f'{symb}'

        return rep

    def _bond_representation(key, just_seen=None):
        key1 = just_seen
        key2 = key
        if key1 is None or key2 is None:
            rep = ''
        else:
            bnd_ord = bnd_ord_dct[frozenset({key1, key2})]
            if bnd_ord == 1:
                rep = ''
            elif bnd_ord == 2:
                rep = '='
            elif bnd_ord == 3:
                rep = '#'
            else:
                raise ValueError("Bond orders greater than 3 not permitted.")
        return rep

    # Get the pool of rings for the graph and set up a dictionary for storing
    # their tags. As the SMILES is built, each next ring that is encountered
    # will be given a tag, removed from the pool, and transferred to the tag
    # dictionary.
    rng_pool = list(rings_atom_keys(rgr))
    rng_tag_dct = {}

    # Determine neighboring keys
    nkeys_dct = dict_.transform_values(atoms_neighbor_atom_keys(rgr), sorted)

    def _recurse_smiles(smi, lst, key, just_seen=None):
        nkeys = nkeys_dct.pop(key) if key in nkeys_dct else []

        # Remove keys just seen from the list of neighbors, to avoid doubling
        # back.
        if just_seen in nkeys:
            nkeys.remove(just_seen)

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
                tags.append(tag)
            if key == rng[0]:
                nkeys.remove(rng[-1])
                closures.append(rng[-1])
                tags.append(tag)

        # Start the SMILES string and connection list. The connection list is
        # used for sorting.
        brep = _bond_representation(key, just_seen)
        arep = _atom_representation(key, just_seen, nkeys, closures=closures)
        rrep = ''.join(map(str, tags))
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

                # # Sort the list of branches by length and index values.
                # srt_idxs = sorted(
                #     range(len(sub_lsts)),
                #     key=lambda i: (len(sub_lsts[i]), sub_lsts[i]))

                # # Apply the sort to both SMILES and lists.
                # sub_smis = list(map(sub_smis.__getitem__, srt_idxs))
                # sub_lsts = list(map(sub_lsts.__getitem__, srt_idxs))

                # Extend the SMILES string
                smi += (''.join(map("({:s})".format, sub_smis[:-1]))
                        + sub_smis[-1])

                # Append the lists of neighboring branches.
                lst.append(sub_lsts)

        return smi, lst

    # If there are terminal atoms, start from the first one
    atm_keys = atom_keys(rgr)
    term_keys = terminal_heavy_atom_keys(gra)
    start_key = min(term_keys) if term_keys else min(atm_keys)

    smi, _ = _recurse_smiles('', [], start_key)

    return smi


def old_smiles(gra, stereo=True, local_stereo=False):
    """ SMILES string from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param stereo: Include stereo?
        :type stereo: bool
        :param local_stereo: Is the graph using local stereo assignments? That
            is, are they based on atom keys rather than canonical keys?
        :type local_stereo: bool
        :returns: the SMILES string
        :rtype: str
    """
    assert is_connected(gra), (
        "Cannot form connection layer for disconnected graph.")

    if not stereo:
        gra = without_stereo_parities(gra)

    # Convert to implicit graph
    gra = implicit(gra)

    # If not using local stereo assignments, canonicalize the graph first.
    # From this point on, the stereo parities can be assumed to correspond to
    # the neighboring atom keys.
    if not local_stereo:
        gra = canonical(gra)

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

    def _atom_representation(key):
        symb = ptab.to_symbol(symb_dct[key])
        if symb in ORGANIC_SUBSET and key not in rad_atm_keys:
            rep = f'{symb}'
        else:
            nhyd = nhyd_dct[key]
            rep = f'[{symb}H{nhyd}]' if nhyd else f'[{symb}]'
        return rep

    def _bond_representation(key1, key2):
        if key1 is None or key2 is None:
            rep = ''
        else:
            bnd_ord = bnd_ord_dct[frozenset({key1, key2})]
            if bnd_ord == 1:
                rep = ''
            elif bnd_ord == 2:
                rep = '='
            elif bnd_ord == 3:
                rep = '#'
            else:
                raise ValueError("Bond orders greater than 3 not permitted.")
        return rep

    # Get the pool of rings for the graph and set up a dictionary for storing
    # their tags. As the SMILES is built, each next ring that is encountered
    # will be given a tag, removed from the pool, and transferred to the tag
    # dictionary.
    rng_pool = list(rings_atom_keys(rgr))
    rng_tag_dct = {}

    # Determine neighboring keys
    nkeys_dct = dict_.transform_values(atoms_neighbor_atom_keys(rgr), list)

    def _recurse_smiles(smi, lst, key, just_seen=None):
        nkeys = nkeys_dct.pop(key) if key in nkeys_dct else []

        # Remove keys just seen from the list of neighbors, to avoid doubling
        # back.
        if just_seen in nkeys:
            nkeys.remove(just_seen)

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
        for rng, tag in rng_tag_dct.items():
            if key == rng[-1]:
                nkeys.remove(rng[0])
                tags.append(tag)
            if key == rng[0]:
                nkeys.remove(rng[-1])
                tags.append(tag)

        # Start the SMILES string and connection list. The connection list is
        # used for sorting.
        brep = _bond_representation(just_seen, key)
        arep = _atom_representation(key)
        rrep = ''.join(map(str, tags))
        smi = f'{brep}{arep}{rrep}'
        lst = [key]

        # Now, extend the layer/list along the neighboring atoms.
        if nkeys:
            # Build sub-strings/lists by recursively calling this function.
            sub_smis = []
            sub_lsts = []
            for nkey in nkeys:
                sub_smi, sub_lst = _recurse_smiles('', [], nkey, just_seen=key)

                sub_smis.append(sub_smi)
                sub_lsts.append(sub_lst)

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
    term_keys = terminal_heavy_atom_keys(gra)
    start_key = min(term_keys) if term_keys else min(atm_keys)

    smi, _ = _recurse_smiles('', [], start_key)

    return smi


if __name__ == '__main__':
    import automol

    # # Rings:
    # ICH = automol.smiles.inchi('C123C(O1)(CCO2)CNC3')
    # ICH = automol.smiles.inchi('C1=C(O[O])NC1=C[CH2]')
    # ICH = automol.smiles.inchi('C1CC1')

    # Tet atoms:
    ICH = automol.smiles.inchi('N[C@](C)(F)C(=O)O')
    ICH = automol.smiles.inchi('N[C@H](C)C(=O)O')
    ICH = automol.smiles.inchi('C[C@H]1CCCCO1')
    GRA = automol.inchi.graph(ICH)
    SMI = smiles(GRA)
    print(SMI)

    SICH = automol.smiles.inchi(SMI)
    print(ICH)
    print(SICH)
    # print(automol.graph.string(GRA, one_indexed=True))
    print(automol.graph.rings_atom_keys(GRA))
    assert SICH == ICH
