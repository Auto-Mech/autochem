"""RDKit interface for working with SMARTS strings."""

from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import ChemicalReaction


def from_smarts(smarts: str) -> ChemicalReaction:
    """Convert SMARTS string to RDKit reaction object.

    :param smarts: SMARTS string
    :return: RDKit reaction object
    """
    return rdChemReactions.ReactionFromSmarts(smarts)


def shape(rxn: ChemicalReaction) -> tuple[int, int]:
    """Get number of reactants and products in reaction.

    :param rxn: RDKit reaction object
    :return: Reaction shape
    """
    return reactant_count(rxn), product_count(rxn)


def reactant_count(rxn: ChemicalReaction) -> int:
    """Get number of reactants in reaction.

    :param rxn: RDKit reaction object
    :return: Number of reactants
    """
    return rxn.GetNumReactantTemplates()


def product_count(rxn: ChemicalReaction) -> int:
    """Get number of products in reaction.

    :param rxn: RDKit reaction object
    :return: Number of products
    """
    return rxn.GetNumProductTemplates()
