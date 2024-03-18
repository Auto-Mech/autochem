""" molecular descriptor libraries

automol module hierarchy
========================

Terminology note: The basic interconvertible molecular types, or "basic types",
are graph, geom, inchi, smiles, and zmat. Types like formula, inchi_key, and
vmat are not considered basic types because they cannot be converted
*back* to the other basic types.

Level 1: No dependencies; no interdependencies

 - const
 - util
 - error
 - mult
 - formula
 - inchi_key
 - vmat
 - prop
 - embed

Level 2: L1 dependencies; hierarchical interdependency (descending)

*The base modules contain functions that do not require convertion to another
basic type.*

 - amchi.base
 - smiles.base
 - geom.base
 - graph.base   [L2 dependencies: geom.base, amchi.base]
 - zmat.base    [L2 dependencies: geom.base]

Level 3: L1-2 dependencies; hierarchical interdependency (descending)

 - extern       [contains RDKit interface needed for working with InChIs]
 - inchi.base   [L3 dependencies: extern]
 - chi.base     [L3 dependencies: extern, inchi.base]

Level 4: L1-3 dependencies; hierarchical interdependency (descending)

*The final modules in level 4 contain all contents from their base modules,
along with additional functions requiring conversion to another basic type.*

 - geom
 - zmat         [L4 dependencies: geom]
 - graph        [L4 dependencies: geom, zmat]
 - amchi        [L4 dependencies: graph, geom]
 - inchi        [L4 dependencies: amchi, graph, geom]
 - chi          [L4 dependencies: amchi, inchi, graph, geom]
 - smiles       [L4 dependencies: graph]

Level 5: L1-4 dependencies; hierarchical interdependency (descending)

 - pot
 - etrans
 - combine
 - reac
 - rotors       [L5 dependencies: reac]
 - symm         [L5 dependencies: reac, rotor]
"""

# L1
from automol import const
from automol import util
from automol import error
from automol import mult
from automol import form
from automol import inchi_key
from automol import vmat
from automol import prop
from automol import embed
# L2
# L3
from automol import extern
# L4
from automol import graph
from automol import geom
from automol import amchi
from automol import inchi
from automol import smiles
from automol import zmat
# L5
from automol import etrans
from automol import combine
from automol import reac
from automol import _deprecated
from automol import data
from automol import symm
# type imports
from automol.const import ReactionClass, ReactionSpin, ReactionInfo


__all__ = [
    # L1
    'const',
    'util',
    'error',
    'mult',
    'form',
    'inchi_key',
    'vmat',
    'prop',
    'embed',
    # L2
    # L3
    'extern',
    # L4
    'graph',
    'geom',
    'amchi',
    'inchi',
    'smiles',
    'smiles',
    'zmat',
    # L5
    'etrans',
    'combine',
    'reac',
    '_deprecated',
    'data',
    'symm',
    # type imports
    'ReactionClass', 'ReactionSpin', 'ReactionInfo',
]
