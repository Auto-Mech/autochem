""" graph functions that depend on resonance structure

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
from automol.graph.base._core import atom_unsaturations
from automol.graph.base._core import bond_unsaturations


if __name__ == '__main__':
    GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 2, None)},
           {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})
    print(atom_unsaturations(GRA))
    print(bond_unsaturations(GRA))
