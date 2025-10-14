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
 - form
 - inchi_key
 - vmat
 - prop
 -  embed

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
# L2
# L3
# L4
# L5
from . import (
    _deprecated,
    amchi,
    combine,
    const,
    data,
    embed,
    error,
    etrans,
    extern,
    form,
    geom,
    graph,
    inchi,
    inchi_key,
    mult,
    prop,
    reac,
    smarts,
    smiles,
    symm,
    util,
    vmat,
    zmat,
)

# type imports
from .const import ReactionClass, ReactionInfo, ReactionSpin

__all__ = [
    # L1
    "const",
    "util",
    "error",
    "mult",
    "form",
    "inchi_key",
    "vmat",
    "prop",
    "embed",
    "smarts",
    # L2
    # L3
    "extern",
    # L4
    "graph",
    "geom",
    "amchi",
    "inchi",
    "smiles",
    "smiles",
    "zmat",
    # L5
    "etrans",
    "combine",
    "reac",
    "_deprecated",
    "data",
    "symm",
    # type imports
    "ReactionClass",
    "ReactionSpin",
    "ReactionInfo",
]
