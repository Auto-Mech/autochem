""" Level 2 RSMILES functions (depend on L1)

The parsing functions apply equally well to SMILES or RSMILES strings, so the
documentation simply refers to SMILES strings.
"""
import string
from typing import List

import numpy
import pyparsing as pp

from phydat import ptab

from ... import util

# Not currently dealing with aromatics
# organic atoms
ORGANIC_SUBSET = ('B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I')
ORGANIC_ATOM = pp.Or(ORGANIC_SUBSET)('symb')

# aromatic atoms
AROMATIC_SUBSET = ('c', 'n', 'o', 'p')
AROMATIC_ATOM = pp.Or(AROMATIC_SUBSET)('symb')

# general atoms
ISO = pp.Opt(pp.Word(pp.nums))
SYMBOL = pp.Combine(
    pp.Char(pp.alphas) + pp.Opt(pp.Char(pp.alphas.lower())))
PARITY = pp.Opt(pp.Or(('@', '@@')))
NHYD = pp.Opt(pp.Combine(pp.Char('H') + pp.Opt(pp.Char(pp.nums))))
CHARGE = pp.Opt(pp.Or((pp.OneOrMore('+') + pp.Opt(pp.nums),
                       pp.OneOrMore('-') + pp.Opt(pp.nums))))
GENERAL_ATOM = (pp.Char('[')('bracket') + ISO('iso') + SYMBOL('symb') +
                PARITY('par') + NHYD('nhyd') + CHARGE('charge') + pp.Char(']'))

# atoms
ATOM = pp.Combine(pp.Or((ORGANIC_ATOM, AROMATIC_ATOM, GENERAL_ATOM)))

# bonds
BOND_STR_2_BOND_ORDER = {'-': 1, '=': 2, '#': 3, None: 1}
BONDS = ('-', '=', '#')
BOND = pp.Opt(pp.Or(BONDS))

DIREC_STR_2_BOOL = {'/': True, '\\': False}
DIRECS = ('/', '\\')
DIREC = pp.Opt(pp.Or(DIRECS))

# ring tags
RING_TAG = pp.Char(pp.nums)
RING_CLOSURE = BOND('bond') + DIREC('direc') + RING_TAG('tag')
RING_CLOSURE_COMBINED = pp.Combine(BOND + DIREC + RING_TAG)
RING_CLOSURES = pp.Opt(pp.OneOrMore(RING_CLOSURE_COMBINED))

# chains
ATOM_ENVIRONMENT = (
    BOND('bond') + DIREC('direc') + ATOM('atom') + RING_CLOSURES('ring_clos'))
ATOM_ENVIRONMENT_COMBINED = pp.Combine(ATOM_ENVIRONMENT)
CHAIN = pp.OneOrMore(ATOM_ENVIRONMENT_COMBINED)
SUBCHAINS = pp.ZeroOrMore(pp.nestedExpr('(', ')', CHAIN))

# smiles
SMILES_PARSER = pp.OneOrMore(CHAIN + SUBCHAINS)


# # conversions
def without_resonance_stereo(smi):
    """ Generate a SMILES string without resonance stereo.

        :param smi: SMILES string
        :type smi: str
        :returns: A SMILES without resonance stereo, which will be readable by
            other codes.
        :rtype: str
    """
    smi = smi.replace('=/', '=')
    smi = smi.replace('=\\', '=')
    return smi


def without_stereo(smi):
    """ Generate a SMILES string without resonance stereo.

        :param smi: SMILES string
        :type smi: str
        :returns: A SMILES without resonance stereo, which will be readable by
            other codes.
        :rtype: str
    """
    smi = smi.replace('/', '')
    smi = smi.replace('\\', '')
    smi = smi.replace('@', '')
    smi = smi.replace('@@', '')
    return smi


def reflect(smi):
    """ If this SMILES has chiral centers, reflect them

        :param smi: SMILES string
        :type smi: str
        :returns: A SMILES string with reflected chiral centers
        :rtype: str
    """
    smi = smi.replace('@@', '&')
    smi = smi.replace('@', '@@')
    smi = smi.replace('&', '@')
    return smi


# # split/join
def split(smi):
    """ Split a SMILES string into connected components

        :param smi: SMILES string
        :type smi: str
        :returns: the component SMILES strings
        :rtype: tuple[str]
    """
    return smi.split('.')


def join(smis):
    """ Join multiple SMILES strings into one multi-component SMILES.

        :param smis: SMILES strings
        :type smis: tuple[str]
        :returns: the multi-component SMILES string
        :rtype: str
    """
    return '.'.join(smis)


def reaction(rsmis, psmis):
    """ Join reactant and product SMILES strings into a reaction SMILES.

        :param rsmis: reactant SMILES strings
        :type rsmis: tuple[str]
        :param psmis: product SMILES strings
        :type psmis: tuple[str]
        :returns: the reaction SMILES string
        :rtype: str
    """
    rsmi = '.'.join(rsmis)
    psmi = '.'.join(psmis)
    rxn_smi = f'{rsmi}>>{psmi}'
    return rxn_smi


def is_reaction(smi: str) -> bool:
    """ Is this a reaction SMILES string?

    :param smi: A reaction SMILES string
    :type smi: str
    :returns: `True` if it is, `False` if it isn't
    :rtype: bool
    """
    return ">>" in smi


def reaction_reagents(smi: str) -> (str, str):
    """ Get the reagents from a reaction SMILES string

    Returns `(None, None)` if it isn't a reaction

    :param smi: A reaction SMILES string
    :type smi: str
    :returns: The reactant and product SMILES strings
    :rtype: (str, str)
    """
    if not is_reaction(smi):
        return (None, None)

    smis = smi.split(">>")
    assert len(smis) == 2, f"Invalid reaction SMILES: {smi}"
    rsmi, psmi = smis
    return rsmi, psmi


def reaction_reactant(smi: str) -> str:
    """ Get the reactant from a reaction SMILES string

    Multiple reactants will be returned as one combined SMILES string
    Returns `None` if it isn't a reaction

    :param smi: A reaction SMILES string
    :type smi: str
    :returns: The reactant SMILES string
    :rtype: str
    """
    rsmi, _ = reaction_reagents(smi)
    return rsmi


def reaction_product(smi: str) -> str:
    """ Get the product from a reaction SMILES string

    Multiple products will be returned as one combined SMILES string
    Returns `None` if it isn't a reaction

    :param smi: A reaction SMILES string
    :type smi: str
    :returns: The product SMILES string
    :rtype: str
    """
    _, psmi = reaction_reagents(smi)
    return psmi


def reaction_reactants(smi: str) -> List[str]:
    """ Get the reactants from a reaction SMILES string

    Multiple reactants will in the order they appear in the string
    Returns `None` if it isn't a reaction

    :param smi: A reaction SMILES string
    :type smi: str
    :returns: The reactant SMILES strings
    :rtype: List[str]
    """
    rsmi = reaction_reactant(smi)
    return None if rsmi is None else rsmi.split('.')


def reaction_products(smi: str) -> List[str]:
    """ Get the products from a reaction SMILES string

    Multiple products will in the order they appear in the string
    Returns `None` if it isn't a reaction

    :param smi: A reaction SMILES string
    :type smi: str
    :returns: The product SMILES strings
    :rtype: List[str]
    """
    rsmi = reaction_product(smi)
    return None if rsmi is None else rsmi.split('.')


def reaction_reactants_and_products(smi: str) -> tuple[list[str], list[str]]:
    """ Get the reactants and products from a reaction SMILES string

    Multiple reactants will be returned as one combined SMILES string
    Returns `None` if it isn't a reaction

    :param smi: A reaction SMILES string
    :returns: The SMILES strings for the reactants and products
    """
    return reaction_reactants(smi), reaction_products(smi)


# # properties
def parse_connected_molecule_properties(smi):
    """ Parse all properties from a SMILES string

        :param smi: SMILES string
        :type smi: str
        :returns: atom symbols by atom key, implicit hydrogens by atom key,
            bond orders by bond key
        :rtype: (dict, dict, dict)
    """
    assert '.' not in smi, f"{smi} is not a connected molecule."

    lst = SMILES_PARSER.parseString(smi).asList()

    # property dictionaries
    symb_dct = {}
    bnd_ord_dct = {}
    nhyd_dct = {}
    atm_par_dct = {}
    bnd_par_dct = {}

    # helper dictionaries
    atm_par_info_dct = {}
    # atom parity info:
    #   {key: (parity, source_key, [hkeys], [ring tags], [nkeys])}
    direc_dct = {}
    # bond direction dictionary:
    #   {(key1, key2): direc}
    rng_clos_dct = {}
    # ring closures dictionary:
    #   {tag: ((key1, key2), bnd_ord, direc)}

    def _recurse_parse(lst, source_key=None):
        # Pop the next atom environment string and parse it
        atm_env_str = lst.pop(0)
        atm_env_dct = ATOM_ENVIRONMENT.parseString(atm_env_str).asDict()

        # Determine the next key
        key = max(symb_dct) + 1 if symb_dct else 0

        # Pop any branches following the current atom
        branch_lsts = []
        while lst and isinstance(lst[0], list):
            branch_lsts.append(lst.pop(0))

        # If there's more to the chain, the rest becomes the last branch
        if lst:
            branch_lsts.append(lst)

        # Read the atomic symbol and update symb_dct
        symb = atm_env_dct['atom']['symb']
        symb_dct[key] = symb

        # Read the bond to the source atom and update bnd_ord_dct
        if source_key is not None:
            bnd_key = frozenset({source_key, key})
            bnd_ord_str = (
                atm_env_dct['bond'] if 'bond' in atm_env_dct else None)
            bnd_ord = BOND_STR_2_BOND_ORDER[bnd_ord_str]
            bnd_ord_dct[bnd_key] = bnd_ord

            if 'direc' in atm_env_dct:
                direc = DIREC_STR_2_BOOL[atm_env_dct['direc']]
                direc_dct[(source_key, key)] = direc

        # Explicit hydrogens become explicit hydrogens in the graph
        hkeys = []
        if 'bracket' in atm_env_dct['atom']:
            if 'nhyd' in atm_env_dct['atom']:
                nhyd_str = atm_env_dct['atom']['nhyd']
                nhyd = _hydrogen_count_from_string(nhyd_str)
                hkey = key
                for _ in range(nhyd):
                    hkey += 1
                    symb_dct[hkey] = 'H'
                    bnd_ord_dct[frozenset({key, hkey})] = 1

                    # Save the hydrogen keys for stereo purposes
                    hkeys.append(hkey)

            # Since we've added explicit hydrogens, set the number of implicit
            # hydrogens to zero
            nhyd_dct[key] = 0

        # Read ring closure bonds and update bnd_ord_dct
        rng_tags = []
        rng_clos_strs = []
        if 'ring_clos' in atm_env_dct:
            rng_clos_strs = atm_env_dct['ring_clos']
            for rng_clos_str in rng_clos_strs:
                bnd_ord, direc, rng_tag = (
                    _ring_closure_from_string(rng_clos_str))
                # If this is a new ring, save the key and the bond order
                if rng_tag not in rng_clos_dct:
                    rng_clos_dct[rng_tag] = (key, bnd_ord, direc)
                # If the ring has been seen before and hasn't been closed,
                # close it now.
                elif isinstance(rng_clos_dct[rng_tag][0], int):
                    clos_key, prev_direc, prev_bnd_ord = rng_clos_dct[rng_tag]
                    bnd_key = frozenset({key, clos_key})
                    # In case only one bond order is specified, iterate over
                    # both and choose whichever one isn't None
                    bnd_ord = next((o for o in (prev_bnd_ord, bnd_ord)
                                    if o is not None), 1)
                    bnd_ord_dct[bnd_key] = bnd_ord

                    # Update the ring closure dictionary
                    rng_clos_dct[rng_tag] = (bnd_key, bnd_ord, direc)

                    # In case only one bond direction is specified, iterate
                    # over both and choose whichever one isn't None
                    # First, flip the previousdirection so that they both
                    # correspond to the same atom ordering
                    prev_direc = (None if prev_direc is None
                                  else not prev_direc)
                    direc = next((d for d in (prev_direc, direc)
                                  if d is not None), None)

                    # Update the direction dictionary
                    if direc is not None:
                        direc_dct[(key, clos_key)] = direc
                # If the ring has been seen before and has been closed, this is
                # a new ring. Edit the ring tag give a unique ID.
                else:
                    for char in string.ascii_uppercase:
                        new_tag = f'{rng_tag}{char}'
                        if new_tag not in rng_clos_dct:
                            break
                    rng_clos_dct[rng_tag] = (key, bnd_ord, direc)

                # Save the tags in order for stereo purposes
                rng_tags.append(rng_tag)

        # Fill out atom parity information for this atom
        if 'par' in atm_env_dct['atom']:
            # Read the parity from the SMILES string
            smi_par = (atm_env_dct['atom']['par'] == '@@')

            # Save information for interpreting the parity
            atm_par_info_dct[key] = (smi_par, source_key, hkeys, rng_tags, [])

        # If the source atom has stereo, add this atom to the list of neighbors
        if source_key in atm_par_info_dct:
            if symb != 'H':
                atm_par_info_dct[source_key][-1].append(key)
            else:
                atm_par_info_dct[source_key][-1].append(-numpy.inf)

        # Iterate over branches and recursively call this function
        for branch_lst in branch_lsts:
            _recurse_parse(branch_lst, source_key=key)

    _recurse_parse(lst)

    # Fill in all explicit hydrogens to simplify stereo
    keys = list(symb_dct.keys())
    symbs = list(symb_dct.values())
    for key, symb in zip(keys, symbs):
        if key not in nhyd_dct:
            nbnds = sum(o for k, o in bnd_ord_dct.items() if key in k)
            nhyd = ptab.valence(symb) - nbnds
            # Remove one hydrogen for aromatics
            if symb == symb.lower():
                nhyd -= 1
            hkey = max(symb_dct.keys())
            for _ in range(nhyd):
                hkey += 1
                symb_dct[hkey] = 'H'
                bnd_ord_dct[frozenset({key, hkey})] = 1

    # Determine local atom parities and fill in atm_par_dct
    for key, vals in atm_par_info_dct.items():
        smi_par, source_key, hkeys, rng_tags, nkeys = vals

        # process source key
        source_keys = [] if source_key is None else [source_key]

        # set hydrogen key to negative infinity
        hkeys = [-numpy.inf for _ in hkeys]

        # process ring tags
        rng_keys = []
        for rng_tag in rng_tags:
            rng_key, = rng_clos_dct[rng_tag][0] - {key}
            rng_keys.append(rng_key)

        skeys = source_keys + hkeys + rng_keys + nkeys
        atm_par = smi_par ^ util.is_odd_permutation(skeys, sorted(skeys))
        atm_par_dct[key] = atm_par

    # Determine local bond parities and fill in bnd_par_dct
    srt_key_dct = {
        k: (k if s != 'H' else -numpy.inf) for k, s in symb_dct.items()}
    bnd_keys = bnd_ord_dct.keys()
    for key1, key2 in bnd_keys:
        nkey1, direc1 = (
            _neighbor_key_and_direction_from_dict(key1, key2, direc_dct))
        nkey2, direc2 = (
            _neighbor_key_and_direction_from_dict(key2, key1, direc_dct))
        if nkey1 is not None and nkey2 is not None:
            smi_par = direc1 ^ direc2

            nkey1s = _neighbor_keys_from_bond_keys(key1, bnd_keys) - {key2}
            nkey2s = _neighbor_keys_from_bond_keys(key2, bnd_keys) - {key1}

            nmax1 = max(nkey1s, key=srt_key_dct.__getitem__)
            nmax2 = max(nkey2s, key=srt_key_dct.__getitem__)

            assert nkey1 in nkey1s, f"{nkey1} not in {nkey1s}"
            assert nkey2 in nkey2s, f"{nkey2} not in {nkey2s}"

            if not (nmax1 == nkey1) ^ (nmax2 == nkey2):
                loc_par = smi_par
            else:
                loc_par = not smi_par

            bnd_par_dct[frozenset({key1, key2})] = loc_par

    return symb_dct, bnd_ord_dct, atm_par_dct, bnd_par_dct


# helpers
def _ring_closure_from_string(rng_clos_str):
    """ Get the bond order and ring tag of a ring closure

        Examples:
         - '=1' would return 2 (order), 1 (tag)
         - '3' would return None (order), 3 (tag)
    """
    parse_dct = RING_CLOSURE.parseString(rng_clos_str).asDict()
    assert 'tag' in parse_dct, f"'tag' not in {parse_dct}"

    if 'bond' in parse_dct:
        bnd_ord = BOND_STR_2_BOND_ORDER[parse_dct['bond']]
    else:
        bnd_ord = None

    if 'direc' in parse_dct:
        direc = DIREC_STR_2_BOOL[parse_dct['direc']]
    else:
        direc = None

    rng_tag = parse_dct['tag']
    return bnd_ord, direc, rng_tag


def _hydrogen_count_from_string(nhyd_str):
    """ Get the hydrogen count from a SMILES hydrogen count string

        Example: 'H' returns 1 and 'H2' returns 2
    """
    if nhyd_str == 'H':
        nhyd = 1
    else:
        nhyd = int(nhyd_str[1:])
    return nhyd


def _neighbor_keys_from_bond_keys(key, bnd_keys):
    """ Determine neighbor keys of an atom from the bond keys
    """
    nkeys = []
    for bnd_key in bnd_keys:
        if key in bnd_key:
            nkey, = bnd_key - {key}
            nkeys.append(nkey)
    return frozenset(nkeys)


def _neighbor_key_and_direction_from_dict(key1, key2, direc_dct):
    r""" Find nkey and its bond direction to key1 in a line-up of the form
         nkey/key1=key2 or nkey\key1=key2, given key1, key2, and the direction
         dictionary.
    """
    keys, direc = next(((ks, d) for ks, d in direc_dct.items()
                        if key1 in ks and key2 not in ks), (None, None))
    if keys is not None:
        idx = keys.index(key1)
        # If key1 is the second element, nkey is the first element.
        # In this case, keep the direction as is.
        if idx == 1:
            nkey = keys[0]
        # If key1 is the first element, nkey is the second element.
        # In this case, flip the direction to correspond to the reverse
        # order
        else:
            nkey = keys[1]
            direc = not direc
    else:
        nkey = None

    return nkey, direc
