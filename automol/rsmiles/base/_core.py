""" Level 2 RSMILES functions (depend on L1)

The parsing functions apply equally well to SMILES or RSMILES strings, so the
documentation simply refers to SMILES strings.
"""

import pyparsing as pp
# from phydat import ptab

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
BONDS = ('-', '=', '#', '/', '\\')
BOND = pp.Opt(pp.Or(BONDS))

# ring tags
RING_TAGS = pp.Opt(pp.OneOrMore(pp.Char(pp.nums)))

# chains
ATOM_ENVIRONMENT = (BOND('bond') + ATOM('atom') + RING_TAGS('ring_tags'))
ATOM_ENVIRONMENT_COMBINED = pp.Combine(ATOM_ENVIRONMENT)
CHAIN = pp.OneOrMore(ATOM_ENVIRONMENT_COMBINED)
SUBCHAINS = pp.ZeroOrMore(pp.nestedExpr('(', ')', content=CHAIN))

# smiles
SMILES_PARSER = pp.OneOrMore(CHAIN + SUBCHAINS)


# # # properties
# def _parse(smi):
#     """ Determine bonds between backbone atoms in a SMILES string
#     """
#     lst = SMILES_PARSER.parseString(smi).asList()
#
#     symb_dct = {}
#     nhyd_dct = {}
#     bnds = {}
#     # atm_par_dct = {}
#     # bnd_par_dct = {}
#
#     entry = lst.pop(0)
#
#     key = 0
#     symb = entry['atom']['symb']
#     symb_dct[key] = symb
#     # nhyd_dct[key] =
#
#     print(ATOM_ENVIRONMENT.parseString(entry).asDict())
#     print(lst)
#
#     entry = lst.pop(0)
#
#     print(ATOM_ENVIRONMENT.parseString(entry).asDict())
#     print(lst)
#
#     entry = lst.pop(0)
#
#     print(ATOM_ENVIRONMENT.parseString(entry).asDict())
#     print(lst)
#
#     entry = lst.pop(0)
#
#     print(ATOM_ENVIRONMENT.parseString(entry).asDict())
#     print(lst)
#
#     entry = lst.pop(0)
#
#     print(ATOM_ENVIRONMENT.parseString(entry).asDict())
#     print(lst)
#
#     entry = lst.pop(0)
#
#     print(ATOM_ENVIRONMENT.parseString(entry).asDict())
#     print(lst)
#
#
# # helpers
# def _implicit_hydrogen_valence(symb, left_bond=None, right_bonds=()):
#     """ Determine the number of implicit hydrogens on an atom
#     """


if __name__ == '__main__':
    SMI = r'FC=CC=C[CH]O'
    # SMI = 'CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O'
    # SMI = r'[C@](O)(Cl)(F)(Br)'
    # _parse(SMI)
