"""reaction-class specific functionality"""

# base reaction class
# tunneling treatments
# submodules
from . import enum, tunnel

# # constructors
# # getters
# # setters
# # others
from ._0core import (
    Reaction,
    apply_zmatrix_conversion,
    class_,
    from_data,
    from_forward_reverse,
    from_string,
    is_radical_radical,
    mapping,
    product_graphs,
    product_mappings,
    product_structures,
    products_graph,
    products_keys,
    reactant_graphs,
    reactant_mappings,
    reactant_structures,
    reactants_graph,
    reactants_keys,
    reagent_mappings,
    relabel,
    remove_atoms,
    reset_conversion_info,
    reverse_without_recalculating,
    set_products_keys,
    set_reactants_keys,
    set_reaction_class,
    set_structures,
    set_ts_graph,
    string,
    structure_type,
    ts_conversion_info,
    ts_graph,
    ts_structure,
    undo_zmatrix_conversion,
    unique,
    update_structures,
    without_dummy_atoms,
    without_stereo,
    without_structures,
)

# stereo-specific reactions
from ._2stereo import expand_stereo, from_old_string, from_string_transitional, reflect

# finders
from ._3find import (
    additions,
    beta_scissions,
    eliminations,
    find,
    hydrogen_abstractions,
    hydrogen_migrations,
    insertions,
    ring_forming_scissions,
    substitutions,
    trivial,
)

# TS geometries
from ._4struc import clean_ts_structure, reverse, with_structures

# conversions
# # constructors from data types
# # converters to various data types
# # additional data types
# # canonicity
from ._5conv import (
    amchis,
    canonical_enantiomer,
    chis,
    display,
    from_amchis,
    from_chis,
    from_geometries,
    from_graphs,
    from_inchis,
    from_smiles,
    from_zmatrices,
    geometries,
    graphs,
    inchis,
    is_canonical_enantiomer,
    reaction_smiles,
    smiles,
    ts_amchi,
    zmatrices,
)

# scan coordinates
from ._6scan import constraint_coordinates, scan_coordinates, scan_values
from ._7scan_deprecated import (
    build_scan_info,
    constraint_coordinate_names,
    scan_coordinate_name,
)

# comp functions
from ._comp import similar_saddle_point_structure

# TS zmatrices
from ._deprecated import zmatrix_coordinate_names

# reaction products
from ._enum import enumerate_reactions, reaction_info_from_string

# species instability transformations
from ._instab import (
    instability_product_graphs,
    instability_product_inchis,
    instability_product_zmas,
    instability_transformation,
)

# phase space theory
from ._pst import pst_cn, pst_kt

__all__ = [
    # base reaction class
    "Reaction",
    # # constructors
    "from_data",
    "from_forward_reverse",
    "from_string",
    "from_old_string",
    "from_string_transitional",
    "string",
    # # getters
    "ts_graph",
    "reactants_keys",
    "products_keys",
    "class_",
    "ts_structure",
    "reactant_structures",
    "product_structures",
    "structure_type",
    "ts_conversion_info",
    # # setters
    "set_ts_graph",
    "set_reactants_keys",
    "set_products_keys",
    "set_reaction_class",
    "set_structures",
    "reset_conversion_info",
    "update_structures",
    # # others
    "reverse_without_recalculating",
    "mapping",
    "reactant_mappings",
    "product_mappings",
    "reagent_mappings",
    "reactant_graphs",
    "product_graphs",
    "reactants_graph",
    "products_graph",
    "relabel",
    "without_stereo",
    "without_structures",
    "apply_zmatrix_conversion",
    "undo_zmatrix_conversion",
    "without_dummy_atoms",
    "remove_atoms",
    "is_radical_radical",
    "unique",
    # stereo-specific reactions
    "expand_stereo",
    "reflect",
    # finders
    "trivial",
    "hydrogen_migrations",
    "beta_scissions",
    "ring_forming_scissions",
    "eliminations",
    "hydrogen_abstractions",
    "additions",
    "insertions",
    "substitutions",
    "find",
    # TS geometries
    "clean_ts_structure",
    "with_structures",
    "reverse",
    # TS zmatrices
    "zmatrix_coordinate_names",
    # conversions
    # # constructors from data types
    "from_graphs",
    "from_amchis",
    "from_inchis",
    "from_chis",
    "from_smiles",
    "from_geometries",
    "from_zmatrices",
    # # converters to various data types
    "graphs",
    "amchis",
    "inchis",
    "chis",
    "smiles",
    "geometries",
    "zmatrices",
    # # additional data types
    "ts_amchi",
    "reaction_smiles",
    "display",
    # # canonicity
    "is_canonical_enantiomer",
    "canonical_enantiomer",
    # scan coordinates
    "scan_coordinates",
    "scan_values",
    "constraint_coordinates",
    "build_scan_info",
    "scan_coordinate_name",
    "constraint_coordinate_names",
    # reaction products
    "enumerate_reactions",
    "reaction_info_from_string",
    # species instability transformations
    "instability_product_zmas",
    "instability_product_inchis",
    "instability_product_graphs",
    "instability_transformation",
    # phase space theory
    "pst_kt",
    "pst_cn",
    # tunneling treatments
    "tunnel",
    # comp functions
    "similar_saddle_point_structure",
    # submodules
    "enum",
]
