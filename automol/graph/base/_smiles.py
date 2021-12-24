""" SMILES (Simplified Molecular Input Line Entry System) functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
from phydat import ptab
from automol.util import dict_
from automol.graph.base._core import atom_keys
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import bond_orders
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import terminal_heavy_atom_keys
from automol.graph.base._core import implicit
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._algo import is_connected
from automol.graph.base._algo import rings_atom_keys
from automol.graph.base._resonance import dominant_resonance
from automol.graph.base._resonance import radical_atom_keys_from_resonance
from automol.graph.base._canon import canonical


ORGANIC_SUBSET = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']


def smiles(gra, stereo=True, can=False):
    """ SMILES string from graph

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param can: Canonicalize the graph? If set to True, the graph will be
            canonicalized first. False by default.
        :type can: bool
        :returns: the SMILES string
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

    # Determine neighboring keys
    nkeys_dct = dict_.transform_values(atoms_neighbor_atom_keys(rgr), list)

    # Find rings for the graph
    rng_atm_keys_lst = list(rings_atom_keys(rgr))
    tagged_rngs_dct = {}

    def _recurse_smiles(smi, lst, key, just_seen=None):
        nkeys = nkeys_dct.pop(key) if key in nkeys_dct else []

        # Remove keys just seen from the list of neighbors, to avoid doubling
        # back.
        if just_seen in nkeys:
            nkeys.remove(just_seen)

        # Start the SMILES string and connection list. The connection list is
        # used for sorting.
        brep = _bond_representation(just_seen, key)
        arep = _atom_representation(key)
        print(f'key, arep, brep: {key} {arep} {brep}')
        smi = f'{brep}{arep}'
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
                print(smi)

                # Extend the list
                lst.extend(sub_lst)
            # If there are multiple neighbors, we joint them as
            #   {arep1}({brep2}{arep2}...)({brep3}{arep3}...){brep4}{arep4}...
            else:
                assert len(sub_lsts) > 1
                print('sub_lsts', sub_lsts)
                print('sub_smis', sub_smis)

                # Sort the list of branches by length and index values.
                srt_idxs = sorted(
                    range(len(sub_lsts)),
                    key=lambda i: (len(sub_lsts[i]), sub_lsts[i]))

                # Apply the sort to both SMILES and lists.
                sub_smis = list(map(sub_smis.__getitem__, srt_idxs))
                sub_lsts = list(map(sub_lsts.__getitem__, srt_idxs))

                # Extend the SMILES string
                smi += (''.join(map("({:s})".format, sub_smis[:-1]))
                        + sub_smis[-1])

                # Append the lists of neighboring branches.
                lst.append(sub_lsts)

        print(f'smi {smi}')
        print(f'lst {lst}')

        return smi, lst

    # If there are terminal atoms, start from the first one
    atm_keys = atom_keys(rgr)
    term_keys = terminal_heavy_atom_keys(gra)
    start_key = min(term_keys) if term_keys else min(atm_keys)

    smi, lst = _recurse_smiles('', [], start_key)
    print(lst)

    return smi


if __name__ == '__main__':
    import automol

    ICH = automol.smiles.inchi('C1=C(O[O])NC1=C[CH2]')
    GRA = automol.inchi.graph(ICH)
    SMI = smiles(GRA)
    print(SMI)

    SICH = automol.smiles.inchi(SMI)
    print(ICH)
    print(SICH)
    assert SICH == ICH
