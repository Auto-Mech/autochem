""" Conversion functions
"""

from collections.abc import Sequence

from .. import amchi, geom, graph
from .. import smiles as smiles_
from .. import zmat as zmat_
from ..extern import rdkit_
from ._0core import (
    Reaction,
    product_graphs,
    product_structures,
    reactant_graphs,
    reactant_structures,
    structure_type,
    ts_graph,
    without_stereo,
)
from ._2stereo import reflect
from ._3find import find
from ._4struc import with_structures


# # constructors from data types
def from_graphs(
    rct_gras: Sequence[object],
    prd_gras: Sequence[object],
    stereo: bool = True,
    struc_typ: str | None = None,
    enant: bool = True,
    strained: bool = False,
) -> tuple[Reaction, ...]:
    """Get reaction objects from graphs.

    :param rct_gras: The reactant graphs
    :param prd_gras: The product graphs
    :param stereo: Include stereoassignments?
    :param struc_typ: Add structures of this type; defaults to "geom"
    :param enant: If expanding stereo, include enantiomers?
    :param strained: If expanding stereo, include strained stereoisomers?
    :returns: A series of reaction objects
    """
    rct_gras = tuple(map(graph.explicit, rct_gras))
    prd_gras = tuple(map(graph.explicit, prd_gras))
    rxns = find(rct_gras, prd_gras, stereo=stereo, enant=enant, strained=strained)

    if struc_typ is not None:
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


def from_amchis(
    rct_achs: Sequence[str],
    prd_achs: Sequence[str],
    stereo: bool = True,
    struc_typ: str | None = None,
    enant: bool = True,
    strained: bool = False,
) -> tuple[Reaction, ...]:
    """Get reaction objects from AMChIs.

    :param rct_achs: The reactant AMChI strings
    :param prd_achs: The product AMChI strings
    :param stereo: Include stereoassignments?
    :param struc_typ: Add structures of this type; defaults to None
    :param enant: If expanding stereo, include enantiomers?
    :param strained: If expanding stereo, include strained stereoisomers?
    :returns: A series of reaction objects
    """
    return from_graphs(
        list(map(amchi.graph, rct_achs)),
        list(map(amchi.graph, prd_achs)),
        stereo=stereo,
        struc_typ=struc_typ,
        enant=enant,
        strained=strained,
    )


def from_inchis(
    rct_ichs: Sequence[str],
    prd_ichs: Sequence[str],
    stereo: bool = True,
    struc_typ: str | None = None,
    enant: bool = True,
    strained: bool = False,
) -> tuple[Reaction, ...]:
    """Get reaction objects from InChIs.

    :param rct_ichs: The reactant InChI strings
    :param prd_ichs: The product InChI strings
    :param stereo: Include stereoassignments?
    :param struc_typ: Add structures of this type; defaults to None
    :param enant: If expanding stereo, include enantiomers?
    :param strained: If expanding stereo, include strained stereoisomers?
    :returns: A series of reaction objects
    """
    return from_amchis(
        rct_ichs,
        prd_ichs,
        stereo=stereo,
        struc_typ=struc_typ,
        enant=enant,
        strained=strained,
    )


def from_chis(
    rct_chis: Sequence[str],
    prd_chis: Sequence[str],
    stereo: bool = True,
    struc_typ: str | None = None,
    enant: bool = True,
    strained: bool = False,
) -> tuple[Reaction, ...]:
    """Get reaction objects from ChIs.

    :param rct_chis: The reactant ChI (InChI or AMChI) strings
    :param prd_chis: The product ChI (InChI or AMChI) strings
    :param stereo: Include stereoassignments?
    :param struc_typ: Add structures of this type; defaults to None
    :param enant: If expanding stereo, include enantiomers?
    :param strained: If expanding stereo, include strained stereoisomers?
    :returns: A series of reaction objects
    """
    return from_amchis(
        rct_chis,
        prd_chis,
        stereo=stereo,
        struc_typ=struc_typ,
        enant=enant,
        strained=strained,
    )


def from_smiles(
    rct_smis: Sequence[str],
    prd_smis: Sequence[str],
    stereo: bool = True,
    struc_typ: str | None = None,
    enant: bool = True,
    strained: bool = False,
) -> tuple[Reaction, ...]:
    """Get reaction objects from SMILES.

    :param rct_smis: The reactant SMILES strings
    :param prd_smis: The product SMILES strings
    :param stereo: Include stereoassignments?
    :param struc_typ: Add structures of this type; defaults to None
    :param enant: If expanding stereo, include enantiomers?
    :param strained: If expanding stereo, include strained stereoisomers?
    :returns: A series of reaction objects
    """
    return from_graphs(
        list(map(smiles_.graph, rct_smis)),
        list(map(smiles_.graph, prd_smis)),
        stereo=stereo,
        struc_typ=struc_typ,
        enant=enant,
        strained=strained,
    )


def from_geometries(
    rct_geos: Sequence[object],
    prd_geos: Sequence[object],
    stereo: bool = True,
    struc_typ: str | None = "geom",
    enant: bool = True,
    strained: bool = False,
) -> tuple[Reaction, ...]:
    """Get reaction objects from geometries.

    :param rct_geos: The reactant geometries
    :param prd_geos: The product geometries
    :param stereo: Include stereoassignments?
    :param struc_typ: Add structures of this type; defaults to "geom"
    :param enant: If expanding stereo, include enantiomers?
    :param strained: If expanding stereo, include strained stereoisomers?
    :returns: A series of reaction objects
    """
    rct_gras = tuple(map(geom.graph, rct_geos))
    prd_gras = tuple(map(geom.graph, prd_geos))
    rxns = find(rct_gras, prd_gras, stereo=stereo, enant=enant, strained=strained)

    if struc_typ is not None:
        rxns = tuple(
            with_structures(r, "geom", rct_strucs=rct_geos, prd_strucs=prd_geos)
            for r in rxns
        )
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


def from_zmatrices(
    rct_zmas: Sequence[object],
    prd_zmas: Sequence[object],
    stereo: bool = True,
    struc_typ: str | None = "zmat",
    enant: bool = True,
    strained: bool = False,
) -> tuple[Reaction, ...]:
    """Get reaction objects from z-matrices.

    :param rct_zmas: The reactant z-matrices
    :param prd_zmas: The product z-matrices
    :param stereo: Include stereoassignments?
    :param struc_typ: Add structures of this type; defaults to "zmat"
    :param enant: If expanding stereo, include enantiomers?
    :param strained: If expanding stereo, include strained stereoisomers?
    :returns: A series of reaction objects
    """
    rct_gras = tuple(map(zmat_.graph, rct_zmas))
    prd_gras = tuple(map(zmat_.graph, prd_zmas))
    rxns = find(rct_gras, prd_gras, stereo=stereo, enant=enant, strained=strained)

    if struc_typ is not None:
        rxns = tuple(
            with_structures(r, "zmat", rct_strucs=rct_zmas, prd_strucs=prd_zmas)
            for r in rxns
        )
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


# # converters to various data types
def graphs(rxn: Reaction, stereo: bool = True, shift_keys: bool = False):
    """Convert the reaction object to graphs.

    :param rxn: the reaction object
    :param stereo: Include stereo? defaults to True
    :type stereo: bool, optional
    :param shift_keys: Shift keys after first reagent, to prevent overlap? default False
    :type shift_keys: bool, optional
    :returns: AMChI strings for the reactants and products
    :rtype: (tuple[str], tuple[str])
    """
    rct_gras = reactant_graphs(rxn, stereo=stereo, shift_keys=shift_keys)
    prd_gras = product_graphs(rxn, stereo=stereo, shift_keys=shift_keys)
    return (rct_gras, prd_gras)


def amchis(rxn: Reaction, stereo: bool = True):
    """Convert the reaction object to AMChIs.

    :param rxn: the reaction object
    :param stereo: Include stereo?
    :type stereo: bool
    :returns: AMChI strings for the reactants and products
    :rtype: (tuple[str], tuple[str])
    """
    rct_gras = reactant_graphs(rxn, stereo=stereo)
    prd_gras = product_graphs(rxn, stereo=stereo)
    rct_chis = tuple(graph.amchi(g, stereo=stereo) for g in rct_gras)
    prd_chis = tuple(graph.amchi(g, stereo=stereo) for g in prd_gras)
    return (rct_chis, prd_chis)


def inchis(rxn: Reaction, stereo: bool = True):
    """Convert the reaction object to InChIs.

    :param rxn: the reaction object
    :param stereo: Include stereo?
    :type stereo: bool
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :returns: InChI strings for the reactants and products
    :rtype: (tuple[str], tuple[str])
    """
    rct_ichs = tuple(graph.inchi(gra, stereo=stereo) for gra in reactant_graphs(rxn))
    prd_ichs = tuple(graph.inchi(gra, stereo=stereo) for gra in product_graphs(rxn))
    return (rct_ichs, prd_ichs)


def chis(rxn: Reaction, stereo: bool = True):
    """Convert the reaction object to ChIs.

    :param rxn: the reaction object
    :param stereo: Include stereo?
    :type stereo: bool
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :returns: ChI strings for the reactants and products
    :rtype: (tuple[str], tuple[str])
    """
    rct_chis = tuple(graph.chi(gra, stereo=stereo) for gra in reactant_graphs(rxn))
    prd_chis = tuple(graph.chi(gra, stereo=stereo) for gra in product_graphs(rxn))
    return (rct_chis, prd_chis)


def smiles(rxn: Reaction, stereo=True, res_stereo=True, exp_singles=False):
    """Convert the reaction object to SMILESs.

    :param rxn: the reaction object
    :param stereo: Include stereo?
    :type stereo: bool
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :param exp_singles: Use explicit '-' for single bonds?
    :type exp_singles: bool
    :returns: SMILES strings for the reactants and products
    :rtype: (tuple[str], tuple[str])
    """
    rct_smis = tuple(
        graph.smiles(gra, stereo=stereo, res_stereo=res_stereo, exp_singles=exp_singles)
        for gra in reactant_graphs(rxn)
    )
    prd_smis = tuple(
        graph.smiles(gra, stereo=stereo, res_stereo=res_stereo, exp_singles=exp_singles)
        for gra in product_graphs(rxn)
    )
    return (rct_smis, prd_smis)


def geometries(rxn: Reaction):
    """Convert the reaction object to geometries.

    :param rxn: the reaction object
    :returns: geometries for the reactants and products
    :rtype: (automol geom data structures, automol geom data structures)
    """
    rct_geos = reactant_structures(rxn) if structure_type(rxn) == "geom" else None
    prd_geos = product_structures(rxn) if structure_type(rxn) == "geom" else None

    if rct_geos is None:
        rct_geos = tuple(map(graph.geometry, reactant_graphs(rxn)))

    if prd_geos is None:
        prd_geos = tuple(map(graph.geometry, product_graphs(rxn)))

    return (rct_geos, prd_geos)


def zmatrices(rxn: Reaction):
    """Convert the reaction object to z-matrices.

    :param rxn: the reaction object
    :returns: z-matrices for the reactants and products
    :rtype: (automol zmat data structures, automol zmat data structures)
    """
    rct_zmas = reactant_structures(rxn) if structure_type(rxn) == "zmat" else None
    prd_zmas = product_structures(rxn) if structure_type(rxn) == "zmat" else None

    if rct_zmas is None:
        rct_zmas = tuple(map(geom.zmatrix, map(graph.geometry, reactant_graphs(rxn))))

    if prd_zmas is None:
        prd_zmas = tuple(map(geom.zmatrix, map(graph.geometry, product_graphs(rxn))))

    return (rct_zmas, prd_zmas)


# # additional data types
def ts_amchi(rxn: Reaction, stereo: bool = True) -> str:
    """Get the AMChI for the reaction TS.

    :param rxn: The reaction object
    :type rxn: Reaction
    :param stereo: Include stereo?
    :type stereo: bool
    :returns: AMChI string for the TS
    :rtype: str
    """
    if not stereo:
        rxn = without_stereo(rxn)

    tsg = ts_graph(rxn)
    ts_chi = graph.amchi(tsg)
    return ts_chi


def reaction_smiles(rxn) -> str:
    """Convert the Reaction object to a reaction SMILES string.

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The reaction SMILES
    :rtype: str
    """
    rct_smis, prd_smis = smiles(rxn)
    rxn_smi = smiles_.reaction(rct_smis, prd_smis)
    return rxn_smi


def display(rxn: Reaction, stereo=True, exp=False, label=False, label_dct=None):
    """Display reaction object to IPython using the RDKit visualizer.

    :param rxn: the reaction object
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param exp: Include explicit hydrogens that aren't needed for stereochemistry?
    :type exp: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    """
    rdkit_.turn_3d_visualization_off()
    return graph.display(
        ts_graph(rxn), stereo=stereo, exp=exp, label=label, label_dct=label_dct
    )


# # canonicity
def is_canonical_enantiomer(srxn: Reaction):
    """Does this reaction have a canonical combination of enantiomers?

    :param srxn: a reaction object with stereo assignments
    :type srxn: Reaction
    :returns: Whether or not the reaction is canonical
    :rtype: bool
    """
    rct_chis, prd_chis = amchis(srxn)
    return amchi.is_canonical_enantiomer_reaction(rct_chis, prd_chis)


def canonical_enantiomer(srxn: Reaction):
    """Convert this reaction into a canonical combination of enantiomers

    :param srxn: a reaction object with stereo assignments
    :type srxn: Reaction
    :returns: Whether or not the reaction is canonical
    :rtype: bool
    """
    if not is_canonical_enantiomer(srxn):
        srxn = reflect(srxn)
    return srxn
