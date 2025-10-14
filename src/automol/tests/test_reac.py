"""Test reac."""

import pytest

from automol import chi as chi_
from automol import geom, graph, reac, smiles, zmat

# Stereo-dependent reaction ID
#   [CH]1[C@H]2CC[C@@H]1O2 => C1=C2CC[C@@H]1O2 + [H]
#         *       ^              *    ^
#
# Notes:
#   1. Without stereochemistry, * and ^ are symmetrically equivalent, so that the
#      C-H bond of either one could be broken to yield the product.
#   2. With stereochemistry, this is not the case. Only one of them can yield the product.
#   3. This was causing a bug in our reaction mapping strategy, as only one of
#      the (apparently equivalent) bonds was enumerated for the beta scission,
#      but it was for the wrong (impossible) stereoisomer.
C5H7O_RGEOS = (
    (
        ("C", (-1.818776, 0.000365, -1.870351)),
        ("C", (-0.728737, 1.849057, 0.019341)),
        ("H", (-3.815204, 0.00202, -2.406224)),
        ("O", (-1.43218, -0.00102, 1.941119)),
        ("C", (2.14661, 1.467934, -0.165888)),
        ("H", (-1.480623, 3.753408, 0.344455)),
        ("C", (-0.7282, -1.849316, 0.017276)),
        ("H", (-1.479758, -3.753888, 0.342001)),
        ("C", (2.147088, -1.467251, -0.165725)),
        ("H", (2.945064, 2.30276, -1.89216)),
        ("H", (3.116856, 2.276812, 1.482364)),
        ("H", (3.11582, -2.275557, 1.483743)),
        ("H", (2.947381, -2.302127, -1.891058)),
    ),
)
C5H7O_PGEOS = (
    (
        ("C", (-1.656404, -1.057905, 1.788813)),
        ("C", (-1.32713, 1.452595, 0.301729)),
        ("H", (-1.047337, -1.533683, 3.704721)),
        ("O", (-1.448161, 0.042876, -2.01624)),
        ("C", (1.639596, 1.769402, 0.474459)),
        ("H", (-2.548583, 3.132631, 0.445664)),
        ("C", (-0.411648, -1.765102, -0.368446)),
        ("C", (2.371036, -1.082399, -0.302428)),
        ("H", (2.362823, 3.238507, -0.810887)),
        ("H", (2.320699, 2.182828, 2.397452)),
        ("H", (3.469656, -2.054244, 1.162259)),
        ("H", (3.335331, -1.208571, -2.134068)),
    ),
    (("H", (0.0, 0.0, 0.0)),),
)


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


@pytest.mark.parametrize(
    "rcts_smi,prds_smi,nexp1,nexp2",
    [
        ("FC=CF.[OH]", "F[CH]C(O)F", 2, 4),
    ],
)
def test__expand_stereo(rcts_smi: str, prds_smi: str, nexp1: int, nexp2: int):
    """Test reac.expand_stereo_for_reaction."""
    print("Testing expand_stereo_for_reaction()")
    print(f"{rcts_smi}>>{prds_smi}")
    rct_smis = rcts_smi.split(".")
    prd_smis = prds_smi.split(".")
    rct_gras0 = tuple(map(smiles.graph, rct_smis))
    prd_gras0 = tuple(map(smiles.graph, prd_smis))
    rxn = reac.find(rct_gras0, prd_gras0, stereo=False)[0]
    srxns = reac.expand_stereo(rxn, enant=False)
    assert len(srxns) == nexp1
    srxns = reac.expand_stereo(rxn, enant=True)
    assert len(srxns) == nexp2


@pytest.mark.parametrize(
    "rcts_smi,prds_smi,enant,count",
    [
        ("F/C=C/F.[OH]", "F[CH][C@H](O)F", False, 1),
        ("CCOCC.[OH]", "C[CH]OCC.O", True, 2),
        ("CCOCC.[OH]", "C[CH]OCC.O", False, 1),
        ("CCO[C@H](C)CC.[OH]", "C[CH]O[C@H](C)CC.O", True, 2),
        ("CCO[C@H](C)CC.[OH]", "C[CH]O[C@H](C)CC.O", False, 2),
    ],
)
def test__expand_stereo_for_reaction(
    rcts_smi: str, prds_smi: str, enant: bool, count: int
):
    """Test reac.expand_stereo_for_reaction."""
    print("Testing expand_stereo_for_reaction()")
    print(f"{rcts_smi}>>{prds_smi}")
    rct_smis = rcts_smi.split(".")
    prd_smis = prds_smi.split(".")
    rct_gras0 = tuple(map(smiles.graph, rct_smis))
    prd_gras0 = tuple(map(smiles.graph, prd_smis))
    rxn = reac.find(rct_gras0, prd_gras0, stereo=False)[0]
    srxns = reac.expand_stereo(rxn, enant=enant, rct_gras=rct_gras0, prd_gras=prd_gras0)
    assert len(srxns) == count

    for srxn in srxns:
        rct_gras1 = reac.reactant_graphs(srxn, shift_keys=False)
        prd_gras1 = reac.product_graphs(srxn, shift_keys=False)
        assert rct_gras1 == rct_gras0
        assert prd_gras1 == prd_gras0


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


@pytest.mark.parametrize(
    "rcts_smi,prds_smi",
    [
        ("CCO.C#[C]", "CC[O].C#C"),
        (r"F\N=[C]/F.[C]#C", r"F\N=C(C#C)/F"),
    ],
)
def test__from_datatypes(rcts_smi: str, prds_smi: str):
    """Test reac.from_<datatype>() functions."""
    print("Testing reac.from_<datatype>() functions")
    print(f"{rcts_smi}>>{prds_smi}")
    rct_smis = rcts_smi.split(".")
    prd_smis = prds_smi.split(".")

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


@pytest.mark.parametrize(
    "fml,rgeos,pgeos,nrxns",
    [
        ("C5H7O", C5H7O_RGEOS, C5H7O_PGEOS, 1),
    ],
)
def test__from_geometries(fml, rgeos, pgeos, nrxns):
    print(f"Testing for {fml}")
    rxns = reac.from_geometries(rct_geos=rgeos, prd_geos=pgeos, stereo=True)
    print(len(rxns))
    assert len(rxns) == nrxns, f"{len(rxns)} != {nrxns}"


@pytest.mark.parametrize(
    "rct_smi,prd_smi",
    [
        # UNIMOLECULAR
        # hydrogen migration
        ("CCCO[O]", "[CH2]CCOO"),
        # hydrogen migration (2TS)
        ("CCC[CH2]", "CC[CH]C"),
        # beta scission (stereo-specific)
        ("F[CH][C@H](O)F", r"F/C=C\F.[OH]"),
        # beta scission (weird z-matrix / linear atom relationships
        ("[CH2]/C=[C]/C#C", "C=C=C=C=[CH].[H]"),
        # beta scission (intramolecular)
        ("C[C]1O[C@H]1COO", r"C/C([O])=C\COO"),
        # ring-forming scission (FIXED)
        ("[CH2]CCCOO", "C1CCCO1.[OH]"),
        (r"[CH2]/C=C\[C@@H](CC)OO", "CC[C@H]1OCC=C1.[OH]"),
        # ring-forming scission with spiro atoms
        ("CC1[C](O1)COO", "CC1C2(O1)CO2.[OH]"),
        ("CC1[C](OC1)CCOO", "CC1C2(OC1)CCO2.[OH]"),
        # elimination
        ("CCCCO[O]", "CCC=C.O[O]"),
        # elimination (HONO)
        ("CCCON(=O)=O", "CCC=O.N(=O)O"),
        # BIMOLECULAR
        # hydrogen abstraction
        ("CCO.[CH3]", "[CH2]CO.C"),
        # hydrogen abstraction (sigma)
        ("CCO.C#[C]", "CC[O].C#C"),
        # hydrogen abstraction (radical radical)
        ("CCC.[H]", "CC[CH2].[HH]"),
        # addition
        ("CC[CH2].[O][O]", "CCCO[O]"),
        # addtition (internal / unimolecular)
        ("CC([O])C=C", "CC(O1)C1[CH2]"),  # not bimolecular
        # addition (H + H => H2)
        ("[H].[H]", "[H][H]"),
        # addition (stereo-specific)
        (r"F/C=C\F.[OH]", "F[CH][C@H](O)F"),
        # addition (stereo-specific with ring)
        ("C1C=C1.[OH]", "C1[CH][C@H]1(O)"),
        # addition (vinyl radical)
        (r"F\N=[C]/F.[C]#C", r"F\N=C(C#C)/F"),
        # addition (vinyl and sigma radicals)
        ("FC=[N].[C]#C", r"F/C=N\C#C"),
        # addition (two vinyl radicals) (FIXED)
        (r"F/C=[C]/[H].[H]/[C]=C/F", r"F/C=C\C=C/F"),
        # addition (case 2)
        ("C=CCCCCCC.[CH2]C", "CCC[CH]CCCCCC"),
        # addition (radical radical 1)
        ("CC[CH2].[H]", "CCC"),
        # addition (radical radical 2) (FIXED)
        ("[H].[OH]", "O"),
        # addition (radical radical 3)
        ("[CH3].[OH]", "CO"),
        # addition (isc??)
        ("N#N.[O]", "[N-]=[N+]=O"),
        # substitution (Sn2) (FIXED)
        ("[C@H](O)(C)F.[Cl]", "[C@@H](O)(C)Cl.[F]"),
        # substitution (FIXED)
        ("CO.[CH2]C", "CCC.[OH]"),
        # insertion
        ("CCC=C.O[O]", "CCCCO[O]"),
        # insertion (HONO)
        ("CCC=O.N(=O)O", "CCCON(=O)=O"),
    ],
)
def test__end_to_end(rct_smi, prd_smi):
    """Test reac.ts_geometry"""
    print("Testing end-to-end functionality")
    print(f"{rct_smi}>>{prd_smi}")

    rct_smis = rct_smi.split(".")
    prd_smis = prd_smi.split(".")

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
    assert reac.without_structures(grxn, keep_info=False) == reac.without_structures(
        grxn_, keep_info=False
    )

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

    # Test that we can build scan information for each class
    scan_names, const_dct, grids, upd_guess = reac.build_scan_info(zrxn, ts_zma)
    print("scan_names:", scan_names)
    print("const_dct:", const_dct)
    print("grids:", grids)
    print("upd_guess:", upd_guess)


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
    # test__from_datatypes("CCO.C#[C]", "CC[O].C#C")
    test__from_datatypes(r"F\N=[C]/F.[C]#C", r"F\N=C(C#C)/F")
    # test__end_to_end("[C@H](O)(C)F.[Cl]", "[C@@H](O)(C)Cl.[F]")
    # test__end_to_end("CCCO[O]", "[CH2]CCOO")
    # test__end_to_end("C[C]1O[C@H]1COO", r"C/C([O])=C\COO")
    # test__end_to_end("CC1[C](O1)COO", "CC1C2(O1)CO2.[OH]")
    # test__end_to_end("CC1[C](O1)COO", "CC1C2(O1)CO2.[OH]")
    # test__end_to_end("CC1[C](OC1)CCOO", "CC1C2(OC1)CCO2.[OH]")
    # test__canonical_enantiomer()
    # test__from_geometries("C5H7O", C5H7O_RGEOS, C5H7O_PGEOS, 1)
    # test__expand_stereo_for_reaction("F/C=C/F.[OH]", "F[CH][C@H](O)F", False, 1)
    # test__expand_stereo_for_reaction("CCOCC.[OH]", "C[CH]OCC.O", True, 2)
    # test__expand_stereo_for_reaction("CCOCC.[OH]", "C[CH]OCC.O", False, 1)
    # test__expand_stereo_for_reaction("CCO[C@H](C)CC.[OH]", "C[CH]O[C@H](C)CC.O", True, 2)
    # test__expand_stereo_for_reaction("CCO[C@H](C)CC.[OH]", "C[CH]O[C@H](C)CC.O", False, 2)
    # test__expand_stereo("FC=CF.[OH]", "F[CH]C(O)F", 2, 4)
