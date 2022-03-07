""" Level 2 RSMILES functions (depend on L1)

The parsing functions apply equally well to SMILES or RSMILES strings, so the
documentation simply refers to SMILES strings.
"""

import pyparsing as pp
from phydat import ptab

# Not currently dealing with aromatics
# organic atoms
ORGANIC_SUBSET = ('B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I')
ORGANIC_ATOM = pp.Or(ORGANIC_SUBSET)('symb')

# general atoms
ISO = pp.Opt(pp.Word(pp.nums))
SYMBOL = pp.Word(pp.alphas, excludeChars='H')
CHIRAL = pp.Opt(pp.Or(('@', '@@')))
NHYD = pp.Opt(pp.Combine('H' + pp.Opt(pp.nums)))
CHARGE = pp.Opt(pp.Or((pp.OneOrMore('+') + pp.Opt(pp.nums),
                       pp.OneOrMore('-') + pp.Opt(pp.nums))))
GENERAL_ATOM = ('[' + ISO('iso') + SYMBOL('symb') + CHIRAL('chi') +
                NHYD('nhyd') + CHARGE('charge') + ']')

# atoms
ATOM = pp.Combine(pp.Or((ORGANIC_ATOM, GENERAL_ATOM)))

# bonds
BOND_STR_2_BOND_ORDER = {'-': 1, '=': 2, '#': 3,
                         None: 1, '/': 1, '\\': 1, '=\\': 2}
BONDS = ('-', '=', '#', '/', '\\', '=\\')
BOND = pp.Opt(pp.Or(BONDS))

# ring tags
RING_TAG = pp.Char(pp.nums)
RING_CLOSURE = BOND('bond') + RING_TAG('tag')
RING_CLOSURE_COMBINED = pp.Combine(BOND + RING_TAG)
RING_CLOSURES = pp.Opt(pp.OneOrMore(RING_CLOSURE_COMBINED))

# chains
ATOM_ENVIRONMENT = (BOND('bond') + ATOM('atom') + RING_CLOSURES('ring_clos'))
ATOM_ENVIRONMENT_COMBINED = pp.Combine(ATOM_ENVIRONMENT)
CHAIN = pp.OneOrMore(ATOM_ENVIRONMENT_COMBINED)
SUBCHAINS = pp.ZeroOrMore(pp.nestedExpr('(', ')', content=CHAIN))

# smiles
SMILES_PARSER = pp.OneOrMore(CHAIN + SUBCHAINS)


# # properties
def parse_properties(smi):
    """ Parse all properties from a SMILES string

        :param smi: SMILES string
        :type smi: str
        :returns: atom symbols by atom key, implicit hydrogens by atom key,
            bond orders by bond key
        :rtype: (dict, dict, dict)
    """
    lst = SMILES_PARSER.parseString(smi).asList()

    # property dictionaries
    symb_dct = {}
    bnd_ord_dct = {}
    nhyd_dct = {}

    # helper dictionaries
    rng_clos_dct = {}  # ring closures: {tag: (bnd_key, bnd_ord)}

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
            bnd_ord = _bond_order_from_string(atm_env_str)
            bnd_ord_dct[bnd_key] = bnd_ord

        # Read ring closure bonds and update bnd_ord_dct
        rng_clos_strs = []
        if 'ring_clos' in atm_env_dct:
            rng_clos_strs = atm_env_dct['ring_clos']
            for rng_clos_str in rng_clos_strs:
                bnd_ord, rng_tag = _ring_closure_from_string(rng_clos_str)
                # If this is a new ring, save the key and the bond order
                if rng_tag not in rng_clos_dct:
                    rng_clos_dct[rng_tag] = (key, bnd_ord)
                # If the ring has been seen before, we are closing it now --
                # save the bond key and bond order. Default to the first
                # explicitly specified bond order.
                else:
                    clos_key, saved_bnd_ord = rng_clos_dct[rng_tag]
                    assert isinstance(clos_key, int), f"{clos_key} not an int"
                    bnd_key = frozenset({key, clos_key})
                    bnd_ord = next((o for o in (saved_bnd_ord, bnd_ord)
                                    if o is not None), 1)
                    bnd_ord_dct[bnd_key] = bnd_ord

        # Explicit hydrogens become explicit hydrogens in the graph
        if 'nhyd' in atm_env_dct['atom']:
            nhyd_str = atm_env_dct['atom']['nhyd']
            nhyd = _hydrogen_count_from_string(nhyd_str)
            hyd_key = key + 1
            for _ in range(nhyd):
                symb_dct[hyd_key] = 'H'
                bnd_ord_dct[frozenset({key, hyd_key})] = 1

            # Since we've added explicit hydrogens, set the number of implicit
            # hydrogens to zero
            nhyd_dct[key] = 0

        # Iterate over branches and recursively call this function
        for branch_lst in branch_lsts:
            _recurse_parse(branch_lst, source_key=key)

    _recurse_parse(lst)

    # Fill in the missing implicit hydrogens
    for key, symb in symb_dct.items():
        if key not in nhyd_dct:
            nbnds = sum(o for k, o in bnd_ord_dct.items() if key in k)
            nhyd_dct[key] = ptab.valence(symb) - nbnds

    return symb_dct, nhyd_dct, bnd_ord_dct


# helpers
def _bond_order_from_string(bnd_str):
    """ Get the bond order from a string that starts with the SMILES encoding
        for bond order

        Example: '=[C@H]...' would return a bond order of 2
    """
    parse_dct = BOND('bond').parseString(bnd_str).asDict()
    bnd_ord_str = parse_dct['bond'] if 'bond' in parse_dct else None
    bnd_ord = BOND_STR_2_BOND_ORDER[bnd_ord_str]
    return bnd_ord


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
    rng_tag = parse_dct['tag']
    return bnd_ord, rng_tag


def _hydrogen_count_from_string(nhyd_str):
    """ Get the hydrogen count from a SMILES hydrogen count string

        Example: 'H' returns 1 and 'H2' returns 2
    """
    if nhyd_str == 'H':
        nhyd = 1
    else:
        nhyd = int(nhyd_str[1:])
    return nhyd


if __name__ == '__main__':
    # SMI = r'FC=CC=C[CH]O'
    # SMI = r'C(=O)C'
    # SMI = r'C=1(OO2)CC2=1'
    SMI = r'C1-C-C=1'
    # SMI = r'C=1(OO2)CC2-1'
    # SMI = r'[C@H](N)(O)(F)'
    # SMI = 'CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O'
    # SMI = r'[C@](O)(Cl)(F)(Br)'
    SYMB_DCT, NHYD_DCT, BND_ORD_DCT = parse_properties(SMI)
    print(SYMB_DCT)
    print(NHYD_DCT)
    print(BND_ORD_DCT)
