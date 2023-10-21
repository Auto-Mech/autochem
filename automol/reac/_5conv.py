""" Conversion functions
"""
from typing import List

import IPython

from automol import amchi, geom, graph
from automol import chi as chi_
from automol import smiles as smiles_
from automol import zmat as zmat_
from automol.extern import rdkit_
from automol.reac._0core import (
    Reaction,
    product_graphs,
    product_structures,
    reactant_graphs,
    reactant_structures,
    structure_type,
    ts_graph,
    without_stereo,
)
from automol.reac._2stereo import reflect
from automol.reac._3find import find
from automol.reac._4struc import with_structures


# # constructors from data types
def from_graphs(
    rct_gras, prd_gras, stereo: bool = True, struc_typ: str = None
) -> List[Reaction]:
    """Get reaction objects from graphs

    :param rct_gras: The reactant graphs
    :type rct_gras: list of automol graphs
    :param prd_gras: The product graphs
    :type prd_gras: list of automol graphs
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to "geom"
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    rct_gras = tuple(map(graph.explicit, rct_gras))
    prd_gras = tuple(map(graph.explicit, prd_gras))
    rxns = find(rct_gras, prd_gras, stereo=stereo)

    if struc_typ is not None:
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


def from_amchis(
    rct_achs, prd_achs, stereo: bool = True, struc_typ: str = None
) -> List[Reaction]:
    """Get reaction objects from AMChIs

    :param rct_achs: The reactant AMChI strings
    :type rct_achs: list[str]
    :param prd_achs: The product AMChI strings
    :type prd_achs: list[str]
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to None
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    return from_chis(rct_achs, prd_achs, stereo=stereo, struc_typ=struc_typ)


def from_inchis(
    rct_ichs, prd_ichs, stereo: bool = True, struc_typ: str = None
) -> List[Reaction]:
    """Get reaction objects from InChIs

    :param rct_ichs: The reactant InChI strings
    :type rct_ichs: list[str]
    :param prd_ichs: The product InChI strings
    :type prd_ichs: list[str]
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to None
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    return from_chis(rct_ichs, prd_ichs, stereo=stereo, struc_typ=struc_typ)


def from_chis(
    rct_chis, prd_chis, stereo: bool = True, struc_typ: str = None
) -> List[Reaction]:
    """Get reaction objects from ChIs

    :param rct_chis: The reactant ChI (InChI or AMChI) strings
    :type rct_chis: list[str]
    :param prd_chis: The product ChI (InChI or AMChI) strings
    :type prd_chis: list[str]
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to None
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    rct_gras = tuple(map(graph.explicit, map(chi_.graph, rct_chis)))
    prd_gras = tuple(map(graph.explicit, map(chi_.graph, prd_chis)))
    rxns = find(rct_gras, prd_gras, stereo=stereo)

    if struc_typ is not None:
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


def from_smiles(
    rct_smis, prd_smis, stereo: bool = True, struc_typ: str = None
) -> List[Reaction]:
    """Get reaction objects from SMILES

    :param rct_smis: The reactant SMILES strings
    :type rct_smis: list[str]
    :param prd_smis: The product SMILES strings
    :type prd_smis: list[str]
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to None
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    rct_gras = tuple(map(graph.explicit, map(smiles_.graph, rct_smis)))
    prd_gras = tuple(map(graph.explicit, map(smiles_.graph, prd_smis)))
    rxns = find(rct_gras, prd_gras, stereo=stereo)

    if struc_typ is not None:
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


def from_geometries(
    rct_geos, prd_geos, stereo: bool = True, struc_typ: str = "geom"
) -> List[Reaction]:
    """Get reaction objects from geometries

    :param rct_geos: The reactant geometries
    :type rct_geos: list of automol geometries
    :param prd_geos: The product geometries
    :type prd_geos: list of automol geometries
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to "geom"
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    rct_gras = tuple(map(geom.graph, rct_geos))
    prd_gras = tuple(map(geom.graph, prd_geos))
    rxns = find(rct_gras, prd_gras, stereo=stereo)

    if struc_typ is not None:
        rxns = tuple(
            with_structures(r, "geom", rct_strucs=rct_geos, prd_strucs=prd_geos)
            for r in rxns
        )
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


def from_zmatrices(
    rct_zmas, prd_zmas, stereo: bool = True, struc_typ: str = "zmat"
) -> List[Reaction]:
    """Get reaction objects from z-matrices

    :param rct_zmas: The reactant z-matrices
    :type rct_zmas: list of automol z-matrices
    :param prd_zmas: The product z-matrices
    :type prd_zmas: list of automol z-matrices
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to "zmat"
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    rct_gras = tuple(map(zmat_.graph, rct_zmas))
    prd_gras = tuple(map(zmat_.graph, prd_zmas))
    rxns = find(rct_gras, prd_gras, stereo=stereo)

    if struc_typ is not None:
        rxns = tuple(
            with_structures(r, "zmat", rct_strucs=rct_zmas, prd_strucs=prd_zmas)
            for r in rxns
        )
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


# # converters to various data types
def graphs(rxn: Reaction, stereo: bool = True, shift_keys: bool = False):
    """Convert the reaction object to graphs

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
    """Convert the reaction object to AMChIs

    :param rxn: the reaction object
    :param stereo: Include stereo?
    :type stereo: bool
    :returns: AMChI strings for the reactants and products
    :rtype: (tuple[str], tuple[str])
    """
    rct_chis = tuple(graph.amchi(gra, stereo=stereo) for gra in reactant_graphs(rxn))
    prd_chis = tuple(graph.amchi(gra, stereo=stereo) for gra in product_graphs(rxn))
    return (rct_chis, prd_chis)


def inchis(rxn: Reaction, stereo: bool = True):
    """Convert the reaction object to InChIs

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
    """Convert the reaction object to ChIs

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
    """Convert the reaction object to SMILESs

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
    """Convert the reaction object to geometries

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
    """Convert the reaction object to z-matrices

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
    """Get the AMChI for the reaction TS

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
    """Convert the Reaction object to a reaction SMILES string

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The reaction SMILES
    :rtype: str
    """
    rct_smis, prd_smis = smiles(rxn)
    rxn_smi = smiles_.reaction(rct_smis, prd_smis)
    return rxn_smi


def rdkit_reaction(rxn: Reaction, stereo=True, res_stereo=False):
    """Convert a reaction object to an RDKit reaction object.

    This is mainly useful for quick visualization with IPython, which can be
    done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_reaction(gra))

        :param rxn: the reaction object
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :param rxn: the reaction object
        :returns: the RDKit reaction
    """
    rdkit_.turn_3d_visualization_off()
    rct_smis, prd_smis = smiles(
        rxn, stereo=stereo, res_stereo=res_stereo, exp_singles=True
    )
    rcts_smi = ".".join(rct_smis)
    prds_smi = ".".join(prd_smis)
    rxn_smi = f"{rcts_smi}>>{prds_smi}"
    return rdkit_.from_smarts(rxn_smi)


def display(rxn: Reaction):
    """Display reaction object to IPython using the RDKit visualizer

    :param rxn: the reaction object
    :returns: None
    """
    rdkit_.turn_3d_visualization_off()
    return IPython.display.display(rdkit_reaction(rxn))


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
