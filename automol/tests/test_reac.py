"""Test automol.reac
"""
import numpy
import automol
from automol import reac


def test__reactant_graphs():
    """Test reac.reactant_graphs"""

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

    _test(["FC=CF", "[OH]"], ["F[CH]C(O)F"])
    _test(["C1CCC1", "[CH3]"], ["C", "C1[CH]CC1"])


def test__expand_stereo():
    """Test reac.expand_stereo_for_reaction"""

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

    _test(["FC=CF", "[OH]"], ["F[CH]C(O)F"], 2, 4)


def test__expand_stereo_for_reaction():
    """Test reac.expand_stereo_for_reaction"""

    def _test(rct_smis, prd_smis):
        print("Testing expand_stereo_for_reaction()")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")
        rct_gras0 = tuple(map(automol.smiles.graph, rct_smis))
        prd_gras0 = tuple(map(automol.smiles.graph, prd_smis))
        rxn = reac.find(rct_gras0, prd_gras0, stereo=False)[0]
        srxns = reac.expand_stereo_for_reaction(rxn, rct_gras0, prd_gras0)
        assert len(srxns) == 1
        (srxn,) = srxns
        rct_gras1 = reac.reactant_graphs(srxn, shift_keys=False)
        prd_gras1 = reac.product_graphs(srxn, shift_keys=False)
        assert rct_gras1 == rct_gras0
        assert prd_gras1 == prd_gras0

    _test(["F/C=C/F", "[OH]"], ["F[CH][C@H](O)F"])


def test__from_old_string():
    """Text reac.from_old_string (with stereo!)"""
    old_rxn_str = """
    reaction class: addition
    forward TS atoms:
        1: {symbol: F, implicit_hydrogen_valence: 0, stereo_parity: null}
        2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
        3: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
        4: {symbol: F, implicit_hydrogen_valence: 0, stereo_parity: null}
        5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
        6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
        7: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
        8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
    forward TS bonds:
        1-2: {order: 1, stereo_parity: null}
        2-3: {order: 1, stereo_parity: true}
        2-5: {order: 1, stereo_parity: null}
        2-7: {order: 0.1, stereo_parity: null}
        3-4: {order: 1, stereo_parity: null}
        3-6: {order: 1, stereo_parity: null}
        7-8: {order: 1, stereo_parity: null}
    reactants keys:
    - [1, 2, 3, 4, 5, 6]
    - [7, 8]
    backward TS atoms:
        1: {symbol: F, implicit_hydrogen_valence: 0, stereo_parity: null}
        2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
        3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
        4: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: false}
        5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
        6: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
        7: {symbol: F, implicit_hydrogen_valence: 0, stereo_parity: null}
        8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
    backward TS bonds:
        1-2: {order: 1, stereo_parity: null}
        2-3: {order: 1, stereo_parity: null}
        2-4: {order: 1, stereo_parity: null}
        4-5: {order: 1, stereo_parity: null}
        4-6: {order: 0.9, stereo_parity: null}
        4-7: {order: 1, stereo_parity: null}
        6-8: {order: 1, stereo_parity: null}
    products keys:
    - [1, 2, 3, 4, 5, 6, 7, 8]
    """
    rxn = reac.from_old_string(old_rxn_str, stereo=True)
    assert rxn == reac.from_data(
        cla="addition",
        rcts_keys=((0, 1, 2, 3, 4, 5), (6, 7)),
        prds_keys=((3, 2, 5, 1, 4, 6, 0, 7),),
        tsg=(
            {
                0: ("F", 0, None),
                1: ("C", 0, False),
                2: ("C", 0, None),
                3: ("F", 0, None),
                4: ("H", 0, None),
                5: ("H", 0, None),
                6: ("O", 0, None),
                7: ("H", 0, None),
            },
            {
                frozenset({0, 1}): (1, None),
                frozenset({1, 2}): (1, True),
                frozenset({1, 4}): (1, None),
                frozenset({1, 6}): (0.1, None),
                frozenset({2, 3}): (1, None),
                frozenset({2, 5}): (1, None),
                frozenset({6, 7}): (1, None),
            },
        ),
    )


def test__find():
    """Test reac.find"""

    def _test(rct_smis, prd_smis):
        print("Testing expand_stereo_for_reaction()")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")
        rct_gras0 = tuple(map(automol.smiles.graph, rct_smis))
        prd_gras0 = tuple(map(automol.smiles.graph, prd_smis))
        srxns = reac.find(rct_gras0, prd_gras0, stereo=True)
        print("len(srxns)", srxns)
        assert len(srxns) > 0
        srxn, *_ = srxns
        rct_gras1 = reac.reactant_graphs(
            srxn, shift_keys=False, standard_reagent_order=False
        )
        prd_gras1 = reac.product_graphs(
            srxn, shift_keys=False, standard_reagent_order=False
        )
        print(f"{rct_gras0}\n---\n{rct_gras1}")
        assert rct_gras1 == rct_gras0
        print(f"{prd_gras0}\n---\n{prd_gras1}")
        assert prd_gras1 == prd_gras0

        # Temporary -- make sure AMChI works for all classes
        ftsg = reac.ts_graph(srxn)
        rtsg = automol.graph.ts.reverse(ftsg)

        fchi = automol.graph.amchi(ftsg)
        rchi = automol.graph.amchi(rtsg)

        print("fchi:", fchi)
        print("rchi:", rchi)
        assert fchi[:-1] == rchi[:-1]

        orig_keys = sorted(automol.graph.atom_keys(ftsg))
        for _ in range(5):
            perm_keys = numpy.random.permutation(orig_keys)
            perm_dct = dict(zip(orig_keys, perm_keys))

            perm_ftsg = automol.graph.relabel(ftsg, perm_dct)
            perm_rtsg = automol.graph.relabel(rtsg, perm_dct)

            perm_fchi = automol.graph.amchi(perm_ftsg)
            perm_rchi = automol.graph.amchi(perm_rtsg)

            print("perm_fchi:", perm_fchi)
            assert perm_fchi == fchi
            print("perm_rchi:", perm_rchi)
            assert perm_rchi == rchi

    # UNIMOLECULAR
    # hydrogen migration
    _test(["CCCO[O]"], ["[CH2]CCOO"])
    # hydrogen migration (2TS)
    _test(["CCC[CH2]"], ["CC[CH]C"])
    # beta scission (stereo-specific)
    _test(["F[CH][C@H](O)F"], [r"F/C=C\F", "[OH]"])
    # ring-forming scission (FIXED)
    _test(["[CH2]CCCOO"], ["C1CCCO1", "[OH]"])
    # elimination
    _test(["CCCCO[O]"], ["CCC=C", "O[O]"])
    # elimination (HONO)
    _test(["CCCON(=O)=O"], ["CCC=O", "N(=O)O"])
    # BIMOLECULAR
    # hydrogen abstraction
    _test(["CCO", "[CH3]"], ["[CH2]CO", "C"])
    # hydrogen abstraction (sigma)
    _test(["CCO", "C#[C]"], ["CC[O]", "C#C"])
    # hydrogen abstraction (radical radical)
    _test(["CCC", "[H]"], ["CC[CH2]", "[HH]"])
    # addition
    _test(["CC[CH2]", "[O][O]"], ["CCCO[O]"])
    # addition (stereo-specific)
    _test([r"F/C=C\F", "[OH]"], ["F[CH][C@H](O)F"])
    # addition (stereo-specific with ring)
    _test(["C1C=C1", "[OH]"], ["C1[CH][C@H]1(O)"])
    # addition (vinyl radical)
    _test([r"F\N=[C]/F", "[C]#C"], [r"F\N=C(C#C)/F"])
    # addition (vinyl and sigma radicals)
    _test(["FC=[N]", "[C]#C"], [r"F/C=N\C#C"])
    # addition (two vinyl radicals) (FIXED)
    _test([r"F/C=[C]/[H]", r"[H]/[C]=C/F"], [r"F/C=C\C=C/F"])
    # addition (case 2)
    _test(["C=CCCCCCC", "[CH2]C"], ["CCC[CH]CCCCCC"])
    # addition (radical radical 1)
    _test(["CC[CH2]", "[H]"], ["CCC"])
    # addition (radical radical 2) (FIXED)
    _test(["[H]", "[OH]"], ["O"])
    # addition (radical radical 3)
    _test(["[CH3]", "[OH]"], ["CO"])
    # addition (isc??)
    _test(["N#N", "[O]"], ["[N-]=[N+]=O"])
    # substitution (Sn2) (FIXED)
    _test(["[C@H](O)(C)F", "[Cl]"], ["[C@@H](O)(C)Cl", "[F]"])
    # substitution (FIXED)
    _test(["CO", "[CH2]C"], ["CCC", "[OH]"])
    # insertion
    _test(["CCC=C", "O[O]"], ["CCCCO[O]"])
    # insertion (HONO)
    _test(["CCC=O", "N(=O)O"], ["CCCON(=O)=O"])


def test__ts_geometry():
    """Test reac.ts_geometry"""

    def _test(rct_smis, prd_smis):
        print("Testing expand_stereo_for_reaction()")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")

        # 1. generate reagent geometries and graphs
        rct_geos = tuple(map(automol.smiles.geometry, rct_smis))
        prd_geos = tuple(map(automol.smiles.geometry, prd_smis))
        rct_gras = tuple(map(automol.geom.graph, rct_geos))
        prd_gras = tuple(map(automol.geom.graph, prd_geos))

        # 2. find reactions
        rxns = reac.find(rct_gras, prd_gras, stereo=True)
        print("len(rxns)", rxns)
        assert len(rxns) > 0
        rxn, *_ = rxns
        rct_gras1 = reac.reactant_graphs(
            rxn, shift_keys=False, standard_reagent_order=False
        )
        prd_gras1 = reac.product_graphs(
            rxn, shift_keys=False, standard_reagent_order=False
        )
        print(f"{rct_gras}\n---\n{rct_gras1}")
        assert rct_gras1 == rct_gras
        print(f"{prd_gras}\n---\n{prd_gras1}")
        assert prd_gras1 == prd_gras

        # 3. generate geometries
        ts_geo = automol.reac.ts_geometry_from_reactants(rxn, rct_geos)
        # Check that the structure is reasonable
        rct_gras_, _ = automol.graph.standard_keys_for_sequence(rct_gras)
        prd_gras_, _ = automol.graph.standard_keys_for_sequence(prd_gras)
        rcts_gra = automol.graph.union_from_sequence(rct_gras_)
        prds_gra = automol.graph.union_from_sequence(prd_gras_)
        rcts_gra = automol.graph.relabel(rcts_gra, automol.reac.mapping(rxn, "R", "T"))
        prds_gra = automol.graph.relabel(prds_gra, automol.reac.mapping(rxn, "P", "T"))
        assert rcts_gra == automol.geom.graph(ts_geo)
        prds_gra1 = automol.graph.set_stereo_from_geometry(prds_gra, ts_geo)
        assert prds_gra == prds_gra1

    # UNIMOLECULAR
    # hydrogen migration
    _test(["CCCO[O]"], ["[CH2]CCOO"])
    # hydrogen migration (2TS)
    _test(["CCC[CH2]"], ["CC[CH]C"])
    # beta scission (stereo-specific)
    _test(["F[CH][C@H](O)F"], [r"F/C=C\F", "[OH]"])
    # ring-forming scission (FIXED)
    _test(["[CH2]CCCOO"], ["C1CCCO1", "[OH]"])
    # elimination
    _test(["CCCCO[O]"], ["CCC=C", "O[O]"])
    # elimination (HONO)
    _test(["CCCON(=O)=O"], ["CCC=O", "N(=O)O"])
    # BIMOLECULAR
    # hydrogen abstraction
    _test(["CCO", "[CH3]"], ["[CH2]CO", "C"])
    # hydrogen abstraction (sigma)
    _test(["CCO", "C#[C]"], ["CC[O]", "C#C"])
    # hydrogen abstraction (radical radical)
    _test(["CCC", "[H]"], ["CC[CH2]", "[HH]"])
    # addition
    _test(["CC[CH2]", "[O][O]"], ["CCCO[O]"])
    # addition (stereo-specific)
    _test([r"F/C=C\F", "[OH]"], ["F[CH][C@H](O)F"])
    # addition (stereo-specific with ring)
    _test(["C1C=C1", "[OH]"], ["C1[CH][C@H]1(O)"])
    # addition (vinyl radical)
    _test([r"F\N=[C]/F", "[C]#C"], [r"F\N=C(C#C)/F"])
    # addition (vinyl and sigma radicals)
    _test(["FC=[N]", "[C]#C"], [r"F/C=N\C#C"])
    # addition (two vinyl radicals) (FIXED)
    _test([r"F/C=[C]/[H]", r"[H]/[C]=C/F"], [r"F/C=C\C=C/F"])
    # addition (case 2)
    _test(["C=CCCCCCC", "[CH2]C"], ["CCC[CH]CCCCCC"])
    # addition (radical radical 1)
    _test(["CC[CH2]", "[H]"], ["CCC"])
    # addition (radical radical 2) (FIXED)
    _test(["[H]", "[OH]"], ["O"])
    # addition (radical radical 3)
    _test(["[CH3]", "[OH]"], ["CO"])
    # addition (isc??)
    _test(["N#N", "[O]"], ["[N-]=[N+]=O"])
    # substitution (Sn2) (FIXED)
    _test(["[C@H](O)(C)F", "[Cl]"], ["[C@@H](O)(C)Cl", "[F]"])
    # substitution (FIXED)
    _test(["CO", "[CH2]C"], ["CCC", "[OH]"])
    # insertion
    _test(["CCC=C", "O[O]"], ["CCCCO[O]"])
    # insertion (HONO)
    _test(["CCC=O", "N(=O)O"], ["CCCON(=O)=O"])


if __name__ == "__main__":
    # test__reactant_graphs()
    # test__expand_stereo()
    # test__expand_stereo_for_reaction()
    # test__from_old_string()
    # test__find()
    test__ts_geometry()
