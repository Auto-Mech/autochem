"""Functions for enumerating reactions."""

from collections.abc import Callable, Mapping, Sequence

from .. import amchi, graph, smiles
from .. import smarts as smarts_
from ._0core import Reaction, from_data


def from_smiles(
    smarts: str, rct_smis: Sequence[str] | Mapping[int, str]
) -> list[Reaction]:
    """Enumerate reactions from SMILES.

    Reactants can be given as lists or dictionaries by position in the SMARTS template.

    :param smarts: Reaction SMARTS string
    :param rct_smis: Reactant SMILES
    :return: List of Reaction objects
    """
    rct_gras = reactants(smarts, rct_smis, conv_=smiles.graph)
    return from_graphs(smarts, rct_gras)


def from_amchis(
    smarts: str, rct_chis: Sequence[str] | Mapping[int, str]
) -> list[Reaction]:
    """Enumerate reactions from AMChIs.

    Reactants can be given as lists or dictionaries by position in the SMARTS template.

    :param rct_chis: Reactant ChIs
    :param smarts: Reaction SMARTS string
    :return: List of Reaction objects
    """
    rct_gras = reactants(smarts, rct_chis, conv_=amchi.graph)
    return from_graphs(smarts, rct_gras)


def from_graphs(
    smarts: str, rct_gras: Sequence[object] | Mapping[int, object]
) -> list[Reaction]:
    """Enumerate reactions from graphs.

    Reactants can be given as lists or dictionaries by position in the SMARTS template.

    :param smarts: SMARTS string
    :param rct_gras: Reactant graphs
    :return: List of Reaction objects
    """
    # Prepare reactants graph
    rct_gras = reactants(smarts, rct_gras, conv_=graph.explicit)
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    rcts_gra = graph.union_from_sequence(rct_gras)

    # Enumerate reactions
    ts_gras = graph.enum.reactions(smarts, rcts_gra)

    # Build reaction objects
    rcts_keys = list(map(sorted, map(graph.atom_keys, rct_gras)))
    rxns = []
    for ts_gra in ts_gras:
        prds_gra = graph.ts.products_graph(ts_gra, stereo=False, dummy=False)
        prd_gras = graph.connected_components(prds_gra)
        prds_keys = list(map(sorted, map(graph.atom_keys, prd_gras)))
        rxns.append(from_data(ts_gra, rcts_keys, prds_keys))
    return rxns


# helpers
def reactants(
    smarts: str,
    rct_vals: Sequence[object] | Mapping[int, object],
    conv_: Callable[[object], object] | None = None,
) -> tuple[list[object], list[object]]:
    """Get lists of reagent values for a SMARTS reaction template.

    :param smarts: SMARTS string
    :param rct_vals: Reagent values
    :return: Lists of reactant and product values
    """
    nrcts = smarts_.reactant_count(smarts)
    rct_vals = dict(enumerate(rct_vals)) if isinstance(rct_vals, Sequence) else rct_vals
    rct_vals = list(map(rct_vals.get, range(nrcts)))
    if conv_ is not None:
        rct_vals = [v if v is None else conv_(v) for v in rct_vals]
    return rct_vals
