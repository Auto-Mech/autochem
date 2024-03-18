"""Test reac
"""
from automol import chi as chi_
from automol import geom, graph, reac, smiles, zmat


def test__reactant_graphs():
    """Test reac.reactant_graphs"""

    def _test(rct_smis, prd_smis):
        print("Testing reactant_graphs()")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")
        rct_gras0 = tuple(map(smiles.graph, rct_smis))
        prd_gras0 = tuple(map(smiles.graph, prd_smis))
        rxns = reac.find(rct_gras0, prd_gras0, stereo=False)
        for rxn in rxns:
            rct_gras1 = reac.reactant_graphs(rxn, shift_keys=False)
            prd_gras1 = reac.product_graphs(rxn, shift_keys=False)
            assert rct_gras1 == rct_gras0
            assert prd_gras1 == prd_gras0

    _test(["FC=CF", "[OH]"], ["F[CH]C(O)F"])
    _test(["C1CCC1", "[CH3]"], ["C", "C1[CH]CC1"])
    _test(["CO", "C[CH2]"], ["CCC", "[OH]"])


def test__expand_stereo():
    """Test reac.expand_stereo_for_reaction"""

    def _test(rct_smis, prd_smis, nexp1, nexp2):
        print("Testing expand_stereo()")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")
        rct_gras0 = tuple(map(smiles.graph, rct_smis))
        prd_gras0 = tuple(map(smiles.graph, prd_smis))
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
        rct_gras0 = tuple(map(smiles.graph, rct_smis))
        prd_gras0 = tuple(map(smiles.graph, prd_smis))
        rxn = reac.find(rct_gras0, prd_gras0, stereo=False)[0]
        srxns = reac.expand_stereo_to_match_reagents(rxn, rct_gras0, prd_gras0)
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


def test__reverse():
    """Test reac.reverse"""

    def _test(rct_smis, prd_smis):
        print("Testing reverse()")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")

        # 1. generate reagent geometries and graphs
        inp_rct_geos = tuple(map(smiles.geometry, rct_smis))
        inp_prd_geos = tuple(map(smiles.geometry, prd_smis))
        inp_rct_gras = tuple(map(geom.graph, inp_rct_geos))
        inp_prd_gras = tuple(map(geom.graph, inp_prd_geos))

        # 2. find reactions
        rxns = reac.find(inp_rct_gras, inp_prd_gras, stereo=True)
        rxn, *_ = rxns  # select the first one for testing

        # 3. add z-matrix structures
        zrxn = reac.with_structures(rxn, "zmat")

        # 4. make sure reversal doesn't break anything
        zrxn0 = reac.reverse(reac.reverse(zrxn))

        # 5. check that we can recover from a `zrxn` without structures, but with dummy
        # atoms
        rxn = reac.without_structures(rxn, keep_info=False)
        rev_rxn = reac.reverse(rxn)
        rev_zrxn = reac.with_structures(rev_rxn, "zmat")
        zrxn1 = reac.reverse(rev_zrxn)

        # 6. tests
        for idx, zrxn_ in enumerate([zrxn0, zrxn1]):
            print(f"Testing z-matrix {idx}")
            ztsg = reac.ts_graph(zrxn_)
            ts_zma = reac.ts_structure(zrxn_)
            rct_zmas = reac.reactant_structures(zrxn_)
            prd_zmas = reac.product_structures(zrxn_)
            rct_zgras = reac.reactant_graphs(zrxn_)
            prd_zgras = reac.product_graphs(zrxn_)

            print(f"\n{ztsg}\n z-matrix matches ? \n{ts_zma}\n")
            assert graph.zmatrix_matches(ztsg, ts_zma)

            assert len(rct_zgras) == len(rct_zmas)
            print("Checking reactant z-matrices....")
            for gra, zma in zip(rct_zgras, rct_zmas):
                print(f"\n{gra}\n z-matrix matches ? \n{zma}\n")
                assert graph.zmatrix_matches(gra, zma)

            assert len(prd_zgras) == len(prd_zmas)
            print("Checking reactant z-matrices....")
            for gra, zma in zip(prd_zgras, prd_zmas):
                print(f"\n{gra}\n z-matrix matches ? \n{zma}\n")
                assert graph.zmatrix_matches(gra, zma)

    _test(["CCO", "C#[C]"], ["CC[O]", "C#C"])


def test__from_datatypes():
    """Test reac.from_<datatype>() functions"""

    def _test(rct_smis, prd_smis):
        print("Testing reac.from_<datatype>() functions")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")

        # 1. generate inputs for various data types
        rct_gras = tuple(map(graph.explicit, map(smiles.graph, rct_smis)))
        prd_gras = tuple(map(graph.explicit, map(smiles.graph, prd_smis)))
        rct_chis = tuple(map(graph.chi, rct_gras))
        prd_chis = tuple(map(graph.chi, prd_gras))
        rct_smis = tuple(map(graph.smiles, rct_gras))
        prd_smis = tuple(map(graph.smiles, prd_gras))
        rct_geos = tuple(map(graph.geometry, rct_gras))
        prd_geos = tuple(map(graph.geometry, prd_gras))
        rct_zmas = tuple(map(geom.zmatrix, rct_geos))
        prd_zmas = tuple(map(geom.zmatrix, prd_geos))

        # 2. get reaction objects from those data types
        rxns_from_gra = reac.from_graphs(rct_gras, prd_gras)
        rxns_from_chi = reac.from_chis(rct_chis, prd_chis)
        rxns_from_smi = reac.from_smiles(rct_smis, prd_smis)
        rxns_from_geo = reac.from_geometries(rct_geos, prd_geos)
        rxns_from_zma = reac.from_zmatrices(rct_zmas, prd_zmas)

        # 3. test the results
        assert (rct_gras, prd_gras) == reac.graphs(rxns_from_gra[0])
        assert (rct_chis, prd_chis) == reac.chis(rxns_from_chi[0])
        assert (rct_smis, prd_smis) == reac.smiles(rxns_from_smi[0])
        assert (rct_geos, prd_geos) == reac.geometries(rxns_from_geo[0])
        assert (rct_zmas, prd_zmas) == reac.zmatrices(rxns_from_zma[0])

    _test(["CCO", "C#[C]"], ["CC[O]", "C#C"])
    _test([r"F\N=[C]/F", "[C]#C"], [r"F\N=C(C#C)/F"])


def test__end_to_end():
    """Test reac.ts_geometry"""

    def _test(rct_smis, prd_smis):
        print("Testing end-to-end functionality")
        print(f"{'.'.join(rct_smis)}>>{'.'.join(prd_smis)}")

        # 1. find reactions
        rxns = reac.from_smiles(rct_smis, prd_smis, stereo=True)
        rxn, *_ = rxns  # select the first one for testing

        # 2. add geometry structures
        grxn = reac.with_structures(rxn, "geom")

        # 3. add z-matrix structures
        zrxn = reac.with_structures(rxn, "zmat")

        # 5. tests
        #   (a.) check that the geometry structures match the reaction graphs
        ts_gra = reac.ts_graph(grxn)
        ts_geo = reac.ts_structure(grxn)
        rct_geos = reac.reactant_structures(grxn)
        prd_geos = reac.product_structures(grxn)
        rct_gras = reac.reactant_graphs(grxn)
        prd_gras = reac.product_graphs(grxn)

        print(f"\n{ts_gra}\n geometry matches ? \n{ts_geo}\n")
        assert graph.geometry_matches(ts_gra, ts_geo)

        assert len(rct_gras) == len(rct_geos)
        print("Checking reactant geometries....")
        for gra, geo in zip(rct_gras, rct_geos):
            print(f"\n{gra}\n geometry matches ? \n{geo}\n")
            assert graph.geometry_matches(gra, geo)

        assert len(prd_gras) == len(prd_geos)
        print("Checking reactant geometries....")
        for gra, geo in zip(prd_gras, prd_geos):
            print(f"\n{gra}\n geometry matches ? \n{geo}\n")
            assert graph.geometry_matches(gra, geo)

        #   (b.) check that the z-matrix structures match the reaction graphs
        ts_zgra = reac.ts_graph(zrxn)
        ts_zma = reac.ts_structure(zrxn)
        rct_zmas = reac.reactant_structures(zrxn)
        prd_zmas = reac.product_structures(zrxn)
        rct_zgras = reac.reactant_graphs(zrxn)
        prd_zgras = reac.product_graphs(zrxn)

        print(f"\n{ts_zgra}\n z-matrix matches ? \n{ts_zma}\n")
        assert graph.zmatrix_matches(ts_zgra, ts_zma)

        assert len(rct_zgras) == len(rct_zmas)
        print("Checking reactant z-matrices....")
        for gra, zma in zip(rct_zgras, rct_zmas):
            print(f"\n{gra}\n z-matrix matches ? \n{zma}\n")
            assert graph.zmatrix_matches(gra, zma)

        assert len(prd_zgras) == len(prd_zmas)
        print("Checking reactant z-matrices....")
        for gra, zma in zip(prd_zgras, prd_zmas):
            print(f"\n{gra}\n z-matrix matches ? \n{zma}\n")
            assert graph.zmatrix_matches(gra, zma)

        #   (c.) check that the z-matrix structure can be converted back to geometries
        grxn_ = reac.with_structures(zrxn, "geom")
        assert reac.without_structures(
            grxn, keep_info=False
        ) == reac.without_structures(grxn_, keep_info=False)

        #   (d.) check that converting to z-matrix again gives the same result
        zrxn_ = reac.with_structures(grxn_, "zmat")
        assert reac.without_structures(zrxn) == reac.without_structures(zrxn_)

        #   (e.) check that we can convert two and from string with structures
        grxn_ = reac.from_string(reac.string(grxn))
        zrxn_ = reac.from_string(reac.string(zrxn))

        print(f"\n{grxn}\n matches ? \n{grxn_}\n")
        assert reac.without_structures(grxn) == reac.without_structures(grxn_)
        assert geom.almost_equal(reac.ts_structure(grxn), reac.ts_structure(grxn_))
        strucs = reac.reactant_structures(grxn) + reac.product_structures(grxn)
        strucs_ = reac.reactant_structures(grxn_) + reac.product_structures(grxn_)
        for struc, struc_ in zip(strucs, strucs_):
            print(f"\n{struc}\n almost equal ? \n{struc_}\n")
            assert geom.almost_equal(struc, struc_)

        print(f"\n{zrxn}\n matches ? \n{zrxn_}\n")
        assert reac.without_structures(zrxn) == reac.without_structures(zrxn_)
        assert zmat.almost_equal(reac.ts_structure(zrxn), reac.ts_structure(zrxn_))
        strucs = reac.reactant_structures(zrxn) + reac.product_structures(zrxn)
        strucs_ = reac.reactant_structures(zrxn_) + reac.product_structures(zrxn_)
        for struc, struc_ in zip(strucs, strucs_):
            print(f"\n{struc}\n almost equal ? \n{struc_}\n")
            assert zmat.almost_equal(struc, struc_)

    # UNIMOLECULAR
    # hydrogen migration
    _test(["CCCO[O]"], ["[CH2]CCOO"])
    # hydrogen migration (2TS)
    _test(["CCC[CH2]"], ["CC[CH]C"])
    # beta scission (stereo-specific)
    _test(["F[CH][C@H](O)F"], [r"F/C=C\F", "[OH]"])
    # ring-forming scission (FIXED)
    _test(["[CH2]CCCOO"], ["C1CCCO1", "[OH]"])
    _test([r"[CH2]/C=C\[C@@H](CC)OO"], ["CC[C@H]1OCC=C1", "[OH]"])
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
    # addtition (internal / unimolecular)
    _test(["CC([O])C=C"], ["CC(O1)C1[CH2]"])  # not bimolecular
    # addition (H + H => H2)
    _test(["[H]", "[H]"], ["[H][H]"])
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


def test__canonical_enantiomer():
    """test reac.canonical_enantiomer"""
    rct_smis = ["CC(OO)C(O[O])C(OO)C"]
    prd_smis = ["CC(OO)C(OO)C(OO)[CH2]"]

    rxns = reac.from_smiles(rct_smis, prd_smis, stereo=False)
    rxn = rxns[0]

    # 2A. Full expansion -- includes non-canonical enantiomer reactions
    print("Full reaction expansion:")
    for srxn in reac.expand_stereo(rxn, enant=True):
        rct_chis, prd_chis = reac.chis(srxn)
        print(" +\n".join(rct_chis) + " =>\n" + " +\n".join(prd_chis))

        # These functions operate directly on the reaction object:
        is_can = reac.is_canonical_enantiomer(srxn)
        print(f"Canonical? {is_can}")
        # Convert it to a canonical enantiomer reaction like this
        srxn = reac.canonical_enantiomer(srxn)
        assert reac.is_canonical_enantiomer(srxn)

        # These are the equivalent functions for ChIs
        is_can = chi_.is_canonical_enantiomer_reaction(rct_chis, prd_chis)
        print(f"Canonical? {is_can}")
        # Convert it to a canonical enantiomer reaction like this
        rct_chis, prd_chis = chi_.canonical_enantiomer_reaction(rct_chis, prd_chis)
        assert chi_.is_canonical_enantiomer_reaction(rct_chis, prd_chis)
        print()

    # 2B. Restricted expansion -- includes only canonical enantiomers
    print("Restricted reaction expansion:")
    for srxn in reac.expand_stereo(rxn, enant=False):
        rct_chis, prd_chis = reac.chis(srxn)
        print(" +\n".join(rct_chis) + " =>\n" + " +\n".join(prd_chis))

        # Check canonicity for a reaction object
        assert reac.is_canonical_enantiomer(srxn)

        # Check canonicity for reaction ChIs
        assert chi_.is_canonical_enantiomer_reaction(rct_chis, prd_chis)
        print()


if __name__ == "__main__":
    # test__reactant_graphs()
    # test__expand_stereo()
    # test__expand_stereo_for_reaction()
    # test__from_old_string()
    # test__reverse()
    # test__from_datatypes()
    test__end_to_end()
    # test__canonical_enantiomer()
