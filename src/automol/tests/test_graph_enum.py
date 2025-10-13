"""test graph.enum."""

import pytest
from automol import graph, smiles


@pytest.mark.parametrize(
    "smarts, rcts_smi, nrxns",
    [(graph.enum.ReactionSmarts.abstraction, "CCC.[OH]", 2)],
)
def test__reactions(smarts, rcts_smi, nrxns):
    """Test reactions."""
    rcts_gra = smiles.graph(rcts_smi)
    ts_gras = graph.enum.reactions(smarts, rcts_gra)
    assert len(ts_gras) == nrxns, f"ts_gras = {ts_gras}"


if __name__ == "__main__":
    test__reactions(graph.enum.ReactionSmarts.abstraction, "CCC.[OH]", 2)
