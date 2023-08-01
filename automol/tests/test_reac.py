"""Test automol.reac
"""
import automol
from automol import reac


def test__reactant_graphs():
    """Test reac.reactant_graphs
    """
    def _test(rct_smis, prd_smis):
        print("Testing reactant_graphs()")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")
        rct_gras0 = tuple(map(automol.smiles.graph, rct_smis))
        prd_gras0 = tuple(map(automol.smiles.graph, prd_smis))
        rxns = reac.find(rct_gras0, prd_gras0, stereo=False)
        for rxn in rxns:
            rct_gras1 = reac.reactant_graphs(rxn, shift_keys=False)
            prd_gras1 = reac.product_graphs(rxn, shift_keys=False)
            assert rct_gras1 == rct_gras0
            assert prd_gras1 == prd_gras0

    _test(['FC=CF', '[OH]'], ['F[CH]C(O)F'])
    _test(['C1CCC1', '[CH3]'], ['C', 'C1[CH]CC1'])


def test__expand_stereo():
    """Test reac.expand_stereo_for_reaction
    """
    def _test(rct_smis, prd_smis, nexp1, nexp2):
        print("Testing expand_stereo()")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")
        rct_gras0 = tuple(map(automol.smiles.graph, rct_smis))
        prd_gras0 = tuple(map(automol.smiles.graph, prd_smis))
        rxn = reac.find(rct_gras0, prd_gras0, stereo=False)[0]
        srxns = reac.expand_stereo(rxn, enant=False)
        assert len(srxns) == nexp1
        srxns = reac.expand_stereo(rxn, enant=True)
        assert len(srxns) == nexp2

    _test(['FC=CF', '[OH]'], ['F[CH]C(O)F'], 2, 4)


def test__expand_stereo_for_reaction():
    """Test reac.expand_stereo_for_reaction
    """
    def _test(rct_smis, prd_smis):
        print("Testing expand_stereo_for_reaction()")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")
        rct_gras0 = tuple(map(automol.smiles.graph, rct_smis))
        prd_gras0 = tuple(map(automol.smiles.graph, prd_smis))
        rxn = reac.find(rct_gras0, prd_gras0, stereo=False)[0]
        srxns = reac.expand_stereo_for_reaction(rxn, rct_gras0, prd_gras0)
        assert len(srxns) == 1
        srxn, = srxns
        rct_gras1 = reac.reactant_graphs(srxn, shift_keys=False)
        prd_gras1 = reac.product_graphs(srxn, shift_keys=False)
        assert rct_gras1 == rct_gras0
        assert prd_gras1 == prd_gras0

    _test(['F/C=C/F', '[OH]'], ['F[CH][C@H](O)F'])


if __name__ == "__main__":
    test__reactant_graphs()
    test__expand_stereo()
    test__expand_stereo_for_reaction()
