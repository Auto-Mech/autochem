""" molecular descriptor libraries

automol module hierarchy
========================

Terminology note: The basic interconvertible molecular types, or "basic types",
are graph, geom, inchi, smiles, and zmat. Types like formula, inchi_key, and
vmat are not considered basic types because they cannot be converted
*back* to the other basic types.

Level 1: No dependencies; no interdependencies

 - par
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

 - geom.base
 - graph.base   [L2 dependencies: geom.base]
 - zmat.base    [L2 dependencies: geom.base]

Level 3: L1-2 dependencies; hierarchical interdependency (descending)

 - extern       [contains RDKit interface needed for working with InChIs]
 - inchi.base   [L3 dependencies: extern]

Level 4: L1-3 dependencies; hierarchical interdependency (descending)

*The final modules in level 4 contain all contents from their base modules,
along with additional functions requiring conversion to another basic type.*

 - geom
 - graph        [L4 dependencies: graph]
 - inchi        [L4 dependencies: graph]
 - smiles
 - zmat         [L4 dependencies: graph, geom]

Level 5: L1-4 dependencies; hierarchical interdependency (descending)

 - pot
 - etrans
 - combine
 - reac
 - rotor        [L5 dependencies: reac]

"""
# L1
from automol import par
from automol import util
from automol import error
from automol import mult
from automol import formula
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
from automol import inchi
from automol import smiles
from automol import zmat
# L5
from automol import pot
from automol import etrans
from automol import combine
from automol import reac
from automol import rotor


__all__ = [
    # L1
    'par',
    'util',
    'error',
    'mult',
    'formula',
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
    'inchi',
    'smiles',
    'zmat',
    # L5
    'pot',
    'etrans',
    'combine',
    'reac',
    'rotor',
]
