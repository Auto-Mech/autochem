""" a molecular graph module

xgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_val1, atm_val2, ...), ...}
bnd_dct: {bnd_key: (bnd_val1, bnd_val2, ...), ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""

# core library
# # constructors
from ..constructors.graph import from_data
from ._core import from_dictionaries
from ._core import from_atoms_and_bonds
from ._core import add_atoms
from ._core import add_bonds
# # value getters
from ._core import atoms
from ._core import bonds
from ._core import atom_keys
from ._core import bond_keys
from ._core import atom_symbols
from ._core import atom_implicit_hydrogen_valences
from ._core import atom_stereo_parities
from ._core import bond_orders
from ._core import bond_stereo_parities
# # value setters
from ._core import relabel
from ._core import set_atom_implicit_hydrogen_valences
from ._core import set_atom_stereo_parities
from ._core import set_bond_orders
from ._core import set_bond_stereo_parities
# # transformations
from ._core import without_bond_orders
from ._core import without_stereo_parities

# graph theory library
# # atom properties
from ._graph import atom_neighbor_keys
from ._graph import atom_bond_keys
from ._graph import atom_neighborhoods
# # bond properties
from ._graph import bond_neighbor_keys
from ._graph import bond_neighborhoods
# # other properties
from ._graph import branch
from ._graph import branch_bond_keys
from ._graph import rings
from ._graph import rings_bond_keys
from ._graph import subgraph
from ._graph import bond_induced_subgraph
# # transformations
from ._graph import delete_atoms

# connectivity graph library
# # atom properties
from ._expl import atom_explicit_hydrogen_valences
from ._expl import atom_explicit_hydrogen_keys
# # other properties
from ._expl import backbone_keys
from ._expl import explicit_hydrogen_keys
# # transformations
from ._expl import add_explicit_hydrogens
from ._expl import implicit
from ._expl import explicit
# # comparisons
from ._expl import backbone_isomorphic
from ._expl import backbone_isomorphism
from ._expl import backbone_unique

# inchi conversion library
from ._inchi import atom_inchi_numbers
from ._inchi import inchi
from ._inchi import stereo_inchi_from_coordinates

# resonance library
# # atom properties
from ._res import atom_bond_valences
from ._res import atom_radical_valences
# # bond properties
from ._res import resonance_dominant_bond_orders
# # other properties
from ._res import maximum_spin_multiplicity
from ._res import possible_spin_multiplicities
# # transformations
from ._res import resonances
from ._res import subresonances
from ._res import dominant_resonances
from ._res import dominant_resonance

# stereo library
# # properties
from ._stereo import stereo_inchi
from ._stereo import is_chiral
from ._stereo import atom_stereo_keys
from ._stereo import bond_stereo_keys
from ._stereo import stereogenic_atom_keys
from ._stereo import stereogenic_bond_keys
# # transformations
from ._stereo import reflection
from ._stereo import stereomers
from ._stereo import substereomers
from ._stereo import enantiomerically_unique

# misc
from ._misc import bond_symmetry_numbers

# submodules
from . import _dict as dict_


__all__ = [
    # core library
    # # constructors
    'from_data',
    'from_dictionaries',
    'from_atoms_and_bonds',
    'add_atoms',
    'add_bonds',
    # # value getters
    'atoms',
    'bonds',
    'atom_keys',
    'bond_keys',
    'atom_symbols',
    'atom_implicit_hydrogen_valences',
    'atom_stereo_parities',
    'bond_orders',
    'bond_stereo_parities',
    # # value setters
    'set_atom_implicit_hydrogen_valences',
    'set_atom_stereo_parities',
    'set_bond_orders',
    'set_bond_stereo_parities',
    # # transformations
    'without_bond_orders',
    'without_stereo_parities',

    # graph theory library
    # # atom properties
    'atom_neighbor_keys',
    'atom_bond_keys',
    'atom_neighborhoods',
    # # bond properties
    'bond_neighbor_keys',
    'bond_neighborhoods',
    # # other properties
    'branch',
    'branch_bond_keys',
    'rings',
    'rings_bond_keys',
    'subgraph',
    'bond_induced_subgraph',
    # # transformations
    'relabel',
    'delete_atoms',

    # chemical connectivity graph library
    # # atom properties
    'atom_explicit_hydrogen_valences',
    'atom_explicit_hydrogen_keys',
    # # other properties
    'backbone_keys',
    'explicit_hydrogen_keys',
    # # transformations
    'add_explicit_hydrogens',
    'implicit',
    'explicit',
    # # comparisons
    'backbone_isomorphic',
    'backbone_isomorphism',
    'backbone_unique',

    # inchi conversion library
    'atom_inchi_numbers',
    'inchi',
    'stereo_inchi_from_coordinates',

    # resonance library
    # # atom properties
    'atom_bond_valences',
    'atom_radical_valences',
    # # bond properties
    'resonance_dominant_bond_orders',
    # # other properties
    'maximum_spin_multiplicity',
    'possible_spin_multiplicities',
    # # transformations
    'resonances',
    'subresonances',
    'dominant_resonances',
    'dominant_resonance',

    # stereo library
    # # properties
    'stereo_inchi',
    'is_chiral',
    'atom_stereo_keys',
    'bond_stereo_keys',
    'stereogenic_atom_keys',
    'stereogenic_bond_keys',
    # # transformations
    'reflection',
    'stereomers',
    'substereomers',
    'enantiomerically_unique',

    # misc
    'bond_symmetry_numbers',

    # submodules
    'dict_',
]
