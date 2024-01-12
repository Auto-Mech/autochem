r""" SMILES (Simplified Molecular Input Line Entry System) functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!


Protocol for encoding resonance single bond stereo into the SMILES string:

    Resonance single bond stereo is specified exactly the same way as double
    bond stereo, using directional bonds on either side of the bond. It isn't
    pretty, but the code is very simple.

    To turn the resonance stereo SMILES into a standard SMILES that is
    interpretable by regular programs, one can simply delete the
    directionality from all double bonds in the string. This may leave some
    extra directional single bonds, but they won't affect anything.
"""
import itertools

import numpy
from phydat import ptab

from automol import util
from automol.graph.base._00core import (
    atom_implicit_hydrogens,
    atom_keys,
    atom_stereo_parities,
    atom_symbols,
    atom_neighbor_atom_key,
    atoms_neighbor_atom_keys,
    bond_orders,
    bond_stereo_parities,
    explicit,
    implicit,
    string,
    terminal_atom_keys,
    with_explicit_stereo_hydrogens,
    without_dummy_atoms,
    without_stereo,
)
from automol.graph.base._02algo import (
    connected_components,
    is_connected,
    rings_atom_keys,
)
from automol.graph.base._03kekule import (
    bad_stereo_bond_keys_from_kekule,
    kekule,
    radical_atom_keys_from_kekule,
)
from automol.graph.base._08canon import canonical
from automol.util import dict_

ORGANIC_SUBSET = ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]

BOND_ORDER_2_BOND_STR = {1: "", 2: "=", 3: "#"}
BOND_ORDER_2_BOND_STR_EXP = {1: "-", 2: "=", 3: "#"}


def smiles(gra, stereo=True, local_stereo=False, res_stereo=True, exp_singles=False):
    """SMILES string from graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereo?
    :type stereo: bool
    :param local_stereo: Is the graph using local stereo assignments? That
        is, are they based on atom keys rather than canonical keys?
    :type local_stereo: bool
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :param exp_singles: Use explicit '-' for single bonds?
    :type exp_singles: bool
    :returns: the SMILES string
    :rtype: str
    """
    gras = connected_components(gra)
    smis = [
        _connected_smiles(
            g,
            stereo=stereo,
            local_stereo=local_stereo,
            res_stereo=res_stereo,
            exp_singles=exp_singles,
        )
        for g in gras
    ]
    smi = ".".join(smis)
    return smi


def _connected_smiles(
    gra, stereo=True, local_stereo=False, res_stereo=True, exp_singles=False
):
    """SMILES string from graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereo?
    :type stereo: bool
    :param local_stereo: Is the graph using local stereo assignments? That
        is, are they based on atom keys rather than canonical keys?
    :type local_stereo: bool
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :param exp_singles: Use explicit '-' for single bonds?
    :type exp_singles: bool
    :returns: the SMILES string
    :rtype: str
    """
    gra = without_dummy_atoms(gra)

    assert is_connected(gra), "Cannot determine SMILES for disconnected graph."

    if not stereo:
        gra = without_stereo(gra)

    # If not using local stereo assignments, canonicalize the graph first.
    # From this point on, the stereo parities can be assumed to correspond to
    # the neighboring atom keys.
    if not local_stereo:
        gra = canonical(gra)

    # Convert to implicit graph
    gra = implicit(gra)

    # Don't allow implicit hydrogens connected to backbone hydrogens
    gra = explicit(gra, atm_keys=atom_keys(gra, symb="H"))

    # Insert hydrogens necessary for bond stereo
    gra = with_explicit_stereo_hydrogens(gra, all_=False, neg=True)

    # Find a dominant resonance
    kgr = kekule(gra, max_stereo_overlap=True)

    # Find stereo parities
    bnd_par_dct = dict_.filter_by_value(
        bond_stereo_parities(kgr), lambda x: x is not None
    )

    # Remove stereo parities if requested
    if not res_stereo:
        bad_bkeys = bad_stereo_bond_keys_from_kekule(kgr)
        bnd_par_dct = dict_.filter_by_key(bnd_par_dct, lambda x: x not in bad_bkeys)

    # Get the pool of stereo bonds for the graph and set up a dictionary for
    # storing bond directions. As the SMILES is built, each stereo
    # bond will be deleted from the pool and the bond directions on either
    # side will be stored. These modifications will be performed by the bond
    # representation function.
    ste_bnd_key_pool = list(bnd_par_dct.keys())
    direc_dct = {}

    # Get the pool of rings for the graph and set up a dictionary for storing
    # their tags. As the SMILES is built, each next ring that is encountered
    # will be given a tag, removed from the pool, and transferred to the tag
    # dictionary. These modifications will be performed by the ring
    # representation function.
    rng_pool = list(rings_atom_keys(kgr))
    rng_tag_dct = {}

    arep_ = atom_representation_generator_(kgr)
    brep_ = bond_representation_generator_(
        kgr, ste_bnd_key_pool, direc_dct, exp_singles=exp_singles
    )
    rrep_ = ring_representation_generator_(
        kgr, direc_dct, rng_pool, rng_tag_dct, exp_singles=exp_singles
    )

    # Determine neighboring keys
    nkeys_dct_pool = dict_.transform_values(atoms_neighbor_atom_keys(kgr), sorted)

    def _recurse_smiles(smi, lst, key, just_seen=None):
        nkeys = nkeys_dct_pool.pop(key) if key in nkeys_dct_pool else []

        # Remove keys just seen from the list of neighbors, to avoid doubling
        # back.
        if just_seen in nkeys:
            nkeys.remove(just_seen)

        # Determine ring, bond, and atom representations for the current atom.
        rrep, nkeys, closures = rrep_(key, nkeys)
        brep = brep_(key, just_seen)
        arep = arep_(key, just_seen, nkeys, closures)

        # Start the SMILES string and connection list. The connection list is
        # used for sorting.
        smi = f"{brep}{arep}{rrep}"
        lst = [key]

        # Now, extend the layer/list along the neighboring atoms.
        if nkeys:
            # Build sub-strings/lists by recursively calling this function.
            sub_smis = []
            sub_lsts = []
            while nkeys:
                nkey = nkeys.pop(0)
                sub_smi, sub_lst = _recurse_smiles("", [], nkey, just_seen=key)

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
                smi += f"{sub_smi}"

                # Extend the list
                lst.extend(sub_lst)
            # If there are multiple neighbors, we joint them as
            #   {arep1}({brep2}{arep2}...)({brep3}{arep3}...){brep4}{arep4}...
            else:
                assert len(sub_lsts) > 1

                # Extend the SMILES string
                smi += "".join(map("({:s})".format, sub_smis[:-1])) + sub_smis[-1]

                # Append the lists of neighboring branches.
                lst.append(sub_lsts)

        return smi, lst

    # If there are terminal atoms, start from the first one
    atm_keys = atom_keys(kgr)
    term_keys = terminal_atom_keys(gra, backbone=False)
    start_key = min(term_keys) if term_keys else min(atm_keys)

    # Don't start on a stereo bond, because this causes confusion of the directionality
    # handling
    ste_keys = set(itertools.chain(*ste_bnd_key_pool))
    if start_key in ste_keys:
        start_key = atom_neighbor_atom_key(gra, start_key, excl_keys=ste_keys)

    smi, _ = _recurse_smiles("", [], start_key)

    return smi


def atom_representation_generator_(kgr):
    """A SMILES atom representation generator.

    SMILES atom representations include the atomic symbol and, if
    applicable, the stereo parity of the atom and the number of hydrogens.

    :param kgr: a kekule graph
    :returns: a function that generates atom representations, with the
        function signature (key, just_seen, nkeys, closures)
    """
    # Determine atom symbols
    symb_dct = atom_symbols(kgr)

    # Determine atom implicit hydrogens
    nhyd_dct = atom_implicit_hydrogens(kgr)

    # Find radical sites for this resonance
    rad_atm_keys = radical_atom_keys_from_kekule(kgr)

    # Find stereo parities
    atm_par_dct = dict_.filter_by_value(
        atom_stereo_parities(kgr), lambda x: x is not None
    )

    def _generator(key, just_seen=None, nkeys=(), closures=()):
        # Determine the atomic symbol.
        symb = ptab.to_symbol(symb_dct[key])

        # Determine the hydrogen count representation.
        nhyd = nhyd_dct[key]
        hyd_rep = f"H{nhyd}" if nhyd > 1 else ("H" if nhyd == 1 else "")

        # Determine the stereo parity representation.
        # If this is not a stereo atom, leave the parity representation blank.
        if key not in atm_par_dct:
            par_rep = ""
        # Otherwise, determine whether the parity is @ (counterclockwise) or @@
        # (clockwise).
        else:
            # Get the list of neighboring keys in order for stereo sorting,
            # including any hydrogens and/or ring closures in the correct
            # order.
            skeys = [just_seen]
            if nhyd:
                assert nhyd == 1
                skeys.append(-numpy.inf)
            if closures:
                skeys.extend(closures)
            skeys.extend(nkeys)

            # Get the canonical parity and compare the SMILES ordering to the
            # canonical ordering to see if the SMILES parity is flipped for
            # this atom (if the permutation is odd).
            can_par = atm_par_dct[key]
            smi_par = can_par ^ util.is_odd_permutation(skeys, sorted(skeys))
            par_rep = "@@" if smi_par else "@"

        # If this is an organic atom with standard valence and no stereo,
        # brackets aren't necessary and the number of hydrogens is implied.
        if (
            symb in ORGANIC_SUBSET
            and key not in rad_atm_keys
            and key not in atm_par_dct
        ):
            rep = f"{symb}"
        # Otherwise, brackets are necessary.
        else:
            rep = f"[{symb}{par_rep}{hyd_rep}]"

        return rep

    return _generator


def bond_representation_generator_(kgr, ste_bnd_key_pool, direc_dct, exp_singles=False):
    r"""A SMILES bond representation generator.

    SMILES bond representations for single, double, and triple bonds are
    given as '', '=', and '#', respectively.

    For stereogenic double bonds, up and down directional single bonds,
    '/' and '\', are used on either side of the double bond to specify its
    stereo parity.

    For stereogenic resonance single bonds, we have our own internal
    system:
        The bond is marked as cis or trans with '^' or '~', respectively.
        The cis/trans specification is with respect to the neighboring
        double bond or directional bond on either side of the resonance
        single bond. If one side has neither a directional bond nor a
        double bond, then the non-directional single bond that is being
        used to specify the stereo parity is marked with a '*'.

    :param kgr: a kekule graph
    :param ste_bnd_key_pool: The pool of stereo bond keys. Each time a new
        one is encountered, it will be removed from the pool.
    :param direc_dct: The dictionary of bond directions. As new stereo
        bonds are encountered, the directions of the bonds on either side
        will be set, and these directions will be saved in this
        dictionary.
    :returns: a function that generates bond representations, with the
        function signature (key, just_seen). The function returns the bond
        representation, along with an updated direc_dct indicating the
        up/down direction that was chosen for the bond, if applicable.
    """
    # Determine bond orders for this resonance
    bnd_ord_dct = bond_orders(kgr)

    # Determine neighbors
    nkeys_dct = atoms_neighbor_atom_keys(kgr)

    # Find stereo parities
    bnd_par_dct = dict_.filter_by_value(
        bond_stereo_parities(kgr), lambda x: x is not None
    )

    def _generator(key, just_seen=None):
        key0 = just_seen
        key1 = key

        # For non-directional bonds, this is all we need:
        if key0 is None:
            # If there is no bond to the previous atom, leave the
            # representation blank.
            rep = ""
        else:
            # If there is a bond to the previous atom, determin the order and
            # set the representation accordingly.
            bnd_ord = bnd_ord_dct[frozenset({key0, key1})]
            rep = (
                BOND_ORDER_2_BOND_STR[bnd_ord]
                if not exp_singles
                else BOND_ORDER_2_BOND_STR_EXP[bnd_ord]
            )

        # Determine if a direction has been assigned to this bond.
        direc = direc_dct[(key0, key1)] if (key0, key1) in direc_dct else ""

        # See if we've encountered a new stereo double bond. This happens when
        # we reach the *first atom* of the stereo double bond.
        ste_bnd_key = next((b for b in ste_bnd_key_pool if key1 in b), None)

        # If we've encountered a new stereo bond, determine the parity, and
        # then determine how directional bonds will be assigned on either
        # side.
        if ste_bnd_key is not None:
            # Remove the new stereo bond from the pool.
            ste_bnd_key_pool.remove(ste_bnd_key)

            # Handle the special case where this is the second directional bond
            # extending from this atom
            same_atom_direc = next(
                (d for k, d in direc_dct.items() if k[1] == just_seen), None
            )
            direc = direc if same_atom_direc is None else same_atom_direc

            # Determine the atoms of the stereo bond and the neighbors that
            # will be assigned directional bonds:
            #   nkey1-key1=key2-nkey2
            (key2,) = ste_bnd_key - {key1}
            nkey1s = set(nkeys_dct[key1]) - {key2}
            nkey2s = set(nkeys_dct[key2]) - {key1}

            nmax1 = max(nkey1s)
            nmax2 = max(nkey2s)

            # The following assumes we have just seen one of the nkey1s, so
            # that the order of bonds is nkey1-key1=key2-nkey2. I think should
            # always be the case.
            assert just_seen in nkey1s, (
                f"{just_seen} not in {nkey1s}. This means we need to "
                f"generalize to the case where there's double bond stereo "
                f"of the form key1(-nkey1)=key2-nkey2."
            )
            # If instead the order were key1(-nkey1)=key2-nkey2, then the
            # flipping rule for cis/trans would be reversed, and we would need
            # to handle things differently because we haven't reached the
            # key1-nkey1 bond yet. I had the code in here to handle both
            # cases, but it was rather hard to parse, so I am trying to keep
            # in simple unless we can't avoid it.

            nkey1 = just_seen
            nkey2 = nmax2

            # Determine parity
            can_par = bnd_par_dct[ste_bnd_key]
            smi_par = can_par if nkey1 == nmax1 else not can_par

            # Set the current bond direction, if it wasn't already set.
            direc = direc if direc else "/"
            flip = not smi_par

            # Set the next bond direction, which should be the same as the
            # first if the bond is trans (parity = True) and should be flipped
            # if the bond is cis (parity = False).
            next_direc = _flip_direction(direc, flip=flip)

            direc_dct[(key2, nkey2)] = next_direc

        if direc and rep == "-":
            rep = direc
        else:
            rep += direc

        return rep

    return _generator


def ring_representation_generator_(
    kgr, direc_dct, rng_pool, rng_tag_dct, exp_singles=False
):
    r"""A SMILES ring representation generator.

    SMILES ring representations for single, double, and triple bonds are
    given as '', '=', and '#', respectively.

    :param kgr: a kekule graph
    :param direc_dct: The dictionary of bond directions used to specify
        bond stereo. This is needed when the directional bond is at the
        end of a ring, so the tag needs to be expressed as /# or \#.
    :param rng_pool: The pool of rings in the graph. Each next ring that
        is encountered will be removed from the pool, given a tag, and
        added to the ring tag dictionary.
    :param rng_tag_dct: The dictionary of assigned numerical tags for
        rings in the graph. These are the numbers used to close the ring,
        for example: C1CCC1 has a ring tag of 1.
    :returns: a function that generates ring representations, with the
        function signature (key, nkeys). The function
        returns the ring representation, along with an updated direc_dct
        indicating the up/down direction that was chosen for the ring, if
        applicable.
    """
    # Determine bond orders for this resonance
    bnd_ord_dct = bond_orders(kgr)

    def _generator(key, nkeys=()):
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
                assert tag < 10, f"Ring tag exceeds 10 for this graph:\n{string(kgr)}"
                rng = util.ring.cycle_item_to_front(new_rng, key, clos_nkey)
                rng_tag_dct[rng] = tag

                # Remove it from the pool of unseen rings
                rng_pool.remove(new_rng)

        tags = []
        closures = []
        for rng, tag in rng_tag_dct.items():
            if key == rng[-1]:
                nkeys.remove(rng[0])
                closures.append(rng[0])
                bnd_ord = bnd_ord_dct[frozenset({rng[-1], rng[0]})]
                bnd_rep = (
                    BOND_ORDER_2_BOND_STR[bnd_ord]
                    if not exp_singles
                    else BOND_ORDER_2_BOND_STR_EXP[bnd_ord]
                )
                # Handle the special case where the last ring bond has stereo
                if (rng[-1], rng[0]) in direc_dct:
                    direc = direc_dct[(rng[-1], rng[0])]
                    bnd_rep = bnd_rep + direc
                tags.append(f"{bnd_rep}{tag}")
            if key == rng[0]:
                nkeys.remove(rng[-1])
                closures.append(rng[-1])
                tags.append(f"{tag}")

        rrep = "".join(map(str, tags))
        return rrep, nkeys, closures

    return _generator


# helpers
def _flip_direction(direc, flip=True):
    """Flip the direction of a directional bond representation

    :param direc: the directional bond representation, '/' or '\\'
    :type direc: str
    :param flip: Flip the representation? If False, don't flip it.
    :type flip: bool
    :returns: the new representation, flipped if requested
    :rtype: str
    """
    if not flip:
        ret = direc
    else:
        ret = "\\" if direc == "/" else "/"
    return ret


if __name__ == "__main__":
    GRA = (
        {
            0: ("H", 0, None),
            1: ("N", 0, None),
            2: ("C", 0, None),
            3: ("C", 0, None),
            4: ("C", 0, None),
            5: ("N", 0, None),
            6: ("H", 0, None),
            7: ("N", 0, None),
            8: ("H", 0, None),
            9: ("H", 0, None),
            10: ("H", 0, None),
        },
        {
            frozenset({3, 4}): (1, None),
            frozenset({10, 4}): (1, None),
            frozenset({2, 3}): (1, None),
            frozenset({3, 7}): (1, False),
            frozenset({1, 2}): (1, False),
            frozenset({4, 5}): (1, True),
            frozenset({0, 1}): (1, None),
            frozenset({8, 7}): (1, None),
            frozenset({5, 6}): (1, None),
            frozenset({9, 2}): (1, None),
        },
    )
    print(smiles(GRA))
