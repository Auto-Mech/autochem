"""Functions for enumerating reactions."""

import itertools
import warnings
from collections.abc import Sequence

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdChemReactions import ChemicalReaction

from ._00core import (
    explicit,
    relabel,
    ts_graph_from_reactants_and_products,
    union_from_sequence,
    without_stereo,
)
from ._02algo import connected_components, unique
from ._03kekule import kekules
from ._12rdkit import from_kekule_graph as to_rdkit
from ._12rdkit import to_graph as from_rdkit

A_ = "Cv4,Ov2"
C_ = "Cv4"
O_ = "Ov2"
As = "CX4,OX2"
Cs = "CX4"
Os = "OX2"
Ar = "Cv3,Ov1"
Cr = "Cv3"
Or = "Ov1"
Au = "CX3,OX1"
Cu = "CX3"
Ou = "OX1"


# Reaction enumeration
class GroupSmarts:
    """SMARTS reaction templates for enumeration."""

    qooh_instability = f"[{Cr}:1][O:2][O:3][H:4]"
    resonant_qooh_instability = f"[*:5]=[*:6][{Cr}:1][O:2][O:3][H:4]"


# Reaction enumeration
class ReactionSmarts:
    """SMARTS reaction templates for enumeration."""

    abstraction = f"[C:1][H:2].[{Or}:3]>>[C:1].[H:2][*:3]"
    h_migration = f"([{Ar}:1].[{A_}:2][H:3])>>([H:3][{A_}:1].[{Ar}:2])"
    h_migration_12 = f"[{Ar}:1][{As}:2][H:3]>>[H:3][{As}:1][{Ar}:2]"
    beta_scission = f"[{A_}:1]-!@[{A_}:2]-[{Ar}:3]>>([{Ar}:1]).([{A_}:2]=[{A_}:3])"
    vinyl_beta_scission = (
        f"[{A_}:1]-!@[{A_}:2]=[{Ar}:3]>>([{Ar}:1]).([{A_}:2]#[{A_}:3])"
    )
    ring_beta_scission = f"[{A_}:1]-@[{A_}:2]-[{Ar}:3]>>([{Ar}:1].[{A_}:2]=[{A_}:3])"
    vinyl_ring_beta_scission = (
        f"[{A_}:1]-@[{A_}:2]=[{Ar}:3]>>([{Ar}:1].[{A_}:2]#[{A_}:3])"
    )
    # Specific
    pi2_addition = f"[*:1]=[*:2].[{Or}:3]>>[*:1]-[*:2]-[*:3]"
    o2_addition = f"[{Ar}:1].[O:2]=[O:3]>>[*:1]-[*:2]-[*:3]"
    qooh_formation = (
        f"([{C_}:1]-[O:2]-[{Or}:3].[{As}:4][H:5])>>([{C_}:1][O:2][O:3][H:5].[{Ar}:4])"
    )
    ho2_elimination = f"[H:5][C:1][C:2][O:3][{Or}:4]>>[C:1]=[C:2].[{Or}:3][O:4][H:5]"
    qooh_beta_scission = (
        f"[{Ar}:1]-[{A_}:2][O:3][O:4][H:5]>>([{Ar}:1]=[{A_}:2]).([O:3][O:4][H:5])"
    )
    qooh_ring_forming_scission = f"([{Cr}:1].[{O_}:2][O:3][H:4])>>[R:1][R:2].[O:3][H:4]"
    qooh_instability = f"[{Cr}:1][O:2][O:3][H:4]>>[C:1]=[O:2].[O:3][H:4]"
    resonant_qooh_instability = (
        f"[*:5]=[*:6][{Cr}:1][O:2][O:3][H:4]>>[*:5]=[*:6][C:1]=[O:2].[O:3][H:4]"
    )


def has_group_match(smarts: str, gra: object) -> bool:
    """Enumerate groups matching a SMARTS string.

    :param smarts: SMARTS pattern for the reaction
    :param gra: Molecular graph representing the reactants
    :return: Graphs
    """
    gra = explicit(gra)
    kgrs = kekules(gra)
    return any(
        itertools.chain.from_iterable(
            group_match_from_kekule(smarts, kgr) for kgr in kgrs
        )
    )


def group_match_from_kekule(smarts: str, kgr: object) -> list[list[int]]:
    """Enumerate matches for a given SMARTS reaction template.

    :param smarts: SMARTS pattern for the reaction
    :param gra: Molecular graph
    :returns: Products graphs
    """
    group = Chem.MolFromSmarts(smarts)
    rdm = to_rdkit(kgr, exp=True, label=True)
    return list(map(list, rdm.GetSubstructMatches(group)))


def reactions(smarts: str, gra: object, symeq: bool = False) -> list[object]:
    """Enumerate reaction TS graphs for a given SMARTS reaction template.

    :param smarts: SMARTS pattern for the reaction
    :param gra: Molecular graph representing the reactants
    :param symeq: Whether to include symmetrically equivalent reactions
    :return: TS graphs
    """
    ts_gras = [
        ts_graph_from_reactants_and_products(gra, p) for p in products(smarts, gra)
    ]
    return ts_gras if symeq else unique(ts_gras)


def products(smarts: str, gra: object) -> list[object]:
    """Enumerate products for a given SMARTS reaction template.

    :param smarts: SMARTS pattern for the reaction
    :param gra: Molecular graph representing the reactants
    :returns: Products graphs
    """
    assert gra == explicit(gra), f"Graph must be explicit\ngra = {gra}"
    gra_without_stereo = without_stereo(gra)
    if gra != gra_without_stereo:
        msg = f"Stereochemistry will be ignored:\n{gra}"
        warnings.warn(msg, stacklevel=2)
        gra = gra_without_stereo

    kgrs = kekules(gra)
    return list(
        itertools.chain.from_iterable(products_from_kekule(smarts, kgr) for kgr in kgrs)
    )


def products_from_kekule(smarts: str, kgr: object) -> list[object]:
    """Enumerate products for a given SMARTS reaction template.

    :param smarts: SMARTS pattern for the reaction
    :param gra: Molecular graph representing the reactants
    :returns: Products graphs
    """
    # Form the reaction object
    rxn = AllChem.ReactionFromSmarts(smarts)

    # Form the reactant graphs
    rct_kgrs = connected_components(kgr)
    # (label=True) stores the current graph keys to the `molAtomMapNumber` property
    rct_rdms = [to_rdkit(kgr, exp=True, label=True) for kgr in rct_kgrs]

    # Enumerate the products
    pos_dct = _template_map_number_to_reactant_position(rxn)
    prds_gras = [
        _products_graph(rct_rdms, prd_rdms, pos_dct)
        for prd_rdms in rxn.RunReactants(rct_rdms)
    ]
    return prds_gras


def _products_graph(
    rct_rdms: Sequence[Mol], prd_rdms: Sequence[Mol], pos_dct: dict[int, int]
) -> object:
    """Map atom keys to reactant positions.

    :param rct_rdms: RDKit reactant molecules
    :param prd_rdms: RDKit product molecules
    :param pos_dct: Mapping of atom map numbers to reactant positions
    :returns: Products graph
    """
    offset_dct = {p + 1: m.GetNumAtoms() for p, m in enumerate(rct_rdms)}

    prd_gras = []
    for prd_pos, prd_rdm in enumerate(prd_rdms):
        key_dct = {}
        for rda in prd_rdm.GetAtoms():
            prop_dct: dict[str, object] = rda.GetPropsAsDict()
            if "molAtomMapNumber" in prop_dct:
                key_dct[rda.GetIdx()] = prop_dct.get("molAtomMapNumber")
            else:
                map_num = rda.GetPropsAsDict().get("old_mapno")
                rct_pos = pos_dct.get(map_num, prd_pos)
                rct_key = rda.GetPropsAsDict().get("react_atom_idx")
                key_dct[rda.GetIdx()] = rct_key + offset_dct.get(rct_pos, 0)
        prd_gra = _product_graph_from_rdkit(prd_rdm)
        prd_gras.append(relabel(prd_gra, key_dct))

    return union_from_sequence(prd_gras)


def _template_map_number_to_reactant_position(rxn: ChemicalReaction) -> dict:
    """Map atom map numbers to reactant positions."""
    map2pos_dct = {}
    for pos in range(rxn.GetNumReactantTemplates()):
        rdt = rxn.GetReactantTemplate(pos)
        for rda in rdt.GetAtoms():
            map2pos_dct[rda.GetAtomMapNum()] = pos
    return map2pos_dct


def _product_graph_from_rdkit(rdm: Mol) -> object:
    """Sanitize an RDKit molecule with radicals."""
    for rda in rdm.GetAtoms():
        rda.SetNoImplicit(True)
    Chem.SanitizeMol(rdm)
    return from_rdkit(rdm)
