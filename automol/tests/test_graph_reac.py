""" test automol.graph
"""

import automol
from automol import graph


def test__trans__is_stereo_compatible():
    """ test graph.trans.is_stereo_compatible
    """
    cgr1 = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
             3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
             6: ('O', 1, None)},
            {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
             frozenset({1, 3}): (1, None)})
    cgr2 = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
             3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
             6: ('O', 1, None)},
            {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({3, 6}): (1, None), frozenset({2, 4}): (1, None),
             frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})

    cgr1 = graph.explicit(cgr1)
    cgr2 = graph.explicit(cgr2)

    ich1 = graph.inchi(cgr1)
    ich2 = graph.inchi(cgr2)
    smi1 = automol.inchi.smiles(ich1)
    smi2 = automol.inchi.smiles(ich2)
    print(smi1)
    print(smi2)

    cgr1s = graph.connected_components(cgr1)
    cgr2s = graph.connected_components(cgr2)
    tras, _, _ = graph.reac.addition(cgr1s, cgr2s)
    for tra in tras:
        assert graph.backbone_isomorphic(graph.trans.apply(tra, cgr1), cgr2)

        sgr1 = graph.stereomers(cgr1)[0]
        print(sgr1)
        for sgr2 in graph.stereomers(cgr2):
            print(sgr2)
            print(graph.trans.is_stereo_compatible(tra, sgr1, sgr2))


def test__reac__hydrogen_migration():
    """ test graph.reac.hydrogen_migration
    """
    # first test a radical site migration
    rct_cgr = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
                3: ('C', 1, None), 4: ('C', 1, None), 5: ('O', 0, None)},
               {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
                frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
                frozenset({0, 1}): (1, None)})

    prd_cgr = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 1, None),
                3: ('C', 1, None), 4: ('C', 0, None), 5: ('O', 0, None)},
               {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
                frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
                frozenset({0, 1}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.hydrogen_migration(rct_cgrs,
                                                             prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    print("hydrogen migration")
    rct_cgr = graph.union_from_sequence(rct_cgrs)
    prd_cgr = graph.union_from_sequence(prd_cgrs)
    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.hydrogen_migration(prd_cgrs,
                                                             rct_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    # then test a tautomerization
    rct_cgr = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('O', 1, None)},
               {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})

    prd_cgr = ({0: ('C', 3, None), 1: ('C', 1, None), 2: ('O', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.hydrogen_migration(rct_cgrs,
                                                             prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.hydrogen_migration(prd_cgrs,
                                                             rct_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs


def test__reac__hydrogen_abstraction():
    """ test graph.reac.hydrogen_abstraction
    """
    rct_cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
                3: ('C', 2, None), 4: ('C', 1, None), 5: ('C', 2, None),
                6: ('C', 2, None), 7: ('O', 1, None)},
               {frozenset({4, 6}): (1, None), frozenset({0, 2}): (1, None),
                frozenset({2, 4}): (1, None), frozenset({5, 6}): (1, None),
                frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})
    prd_cgr = ({0: ('C', 2, None), 1: ('C', 3, None), 2: ('C', 1, None),
                3: ('C', 2, None), 4: ('C', 1, None), 5: ('C', 2, None),
                6: ('C', 2, None), 7: ('O', 2, None)},
               {frozenset({4, 6}): (1, None), frozenset({0, 2}): (1, None),
                frozenset({2, 4}): (1, None), frozenset({5, 6}): (1, None),
                frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.hydrogen_abstraction(rct_cgrs,
                                                               prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    print("hydrogen abstraction")
    rct_cgr = graph.union_from_sequence(rct_cgrs)
    print(rct_cgr)
    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)

    tras, prd_idxs, rct_idxs = graph.reac.hydrogen_abstraction(prd_cgrs,
                                                               rct_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, prd_cgr), rct_cgr)


def test__reac__addition():
    """ test graph.reac.addition
    """
    rct_cgr = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
                3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
                6: ('O', 1, None)},
               {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
                frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
                frozenset({1, 3}): (1, None)})
    prd_cgr = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
                3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
                6: ('O', 1, None)},
               {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
                frozenset({3, 6}): (1, None), frozenset({2, 4}): (1, None),
                frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.addition(rct_cgrs, prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    print("addition")
    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)


def test__reac__beta_scission():
    """ test graph.reac.beta_scission
    """
    rct_cgr = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
                3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
                6: ('O', 1, None)},
               {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
                frozenset({3, 6}): (1, None), frozenset({2, 4}): (1, None),
                frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})
    prd_cgr = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
                3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
                6: ('O', 1, None)},
               {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
                frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
                frozenset({1, 3}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.beta_scission(rct_cgrs, prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    print("beta scission")
    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)


def test__reac__elimination():
    """ test graph.reac.elimination
    """
    # CH3CH3OO => CH2CH2 + HOO
    rct_cgr = ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('O', 0, None),
                3: ('O', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({1, 3}): (1, None),
                frozenset({2, 3}): (1, None)})
    prd_cgr = ({0: ('C', 2, None), 1: ('C', 2, None), 2: ('O', 1, None),
                3: ('O', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.elimination(rct_cgrs, prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    print("elimination")
    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)

    # CH3C(=O)CH3 => CH3CH3 + :C=O
    rct_cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 0, None),
                3: ('O', 0, None)},
               {frozenset({0, 2}): (1, None), frozenset({2, 3}): (1, None),
                frozenset({1, 2}): (1, None)})
    prd_cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 0, None),
                3: ('O', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.elimination(rct_cgrs, prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)

    # CH3C(=O)OCH3 => CH3CH3 + O=C=O
    rct_cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 0, None),
                3: ('O', 0, None), 4: ('O', 0, None)},
               {frozenset({0, 2}): (1, None), frozenset({1, 4}): (1, None),
                frozenset({2, 3}): (1, None), frozenset({2, 4}): (1, None)})
    prd_cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 0, None),
                3: ('O', 0, None), 4: ('O', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
                frozenset({2, 4}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.elimination(rct_cgrs, prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)


def test__reac__insertion():
    """ test graph.reac.insertion
    """
    # CH2CH2 + HOO => CH3CH3OO
    rct_cgr = ({0: ('C', 2, None), 1: ('C', 2, None), 2: ('O', 1, None),
                3: ('O', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None)})
    prd_cgr = ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('O', 0, None),
                3: ('O', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({1, 3}): (1, None),
                frozenset({2, 3}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.insertion(rct_cgrs, prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    print("insertion")
    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)

    # CH3CH3 + :C=O => CH3C(=O)CH3 (CO insertion)
    rct_cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 0, None),
                3: ('O', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None)})
    prd_cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 0, None),
                3: ('O', 0, None)},
               {frozenset({0, 2}): (1, None), frozenset({2, 3}): (1, None),
                frozenset({1, 2}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.insertion(rct_cgrs, prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)

    # CH3CH3 + O=C=O => CH3C(=O)OCH3 (CO2 insertion)
    rct_cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 0, None),
                3: ('O', 0, None), 4: ('O', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
                frozenset({2, 4}): (1, None)})
    prd_cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 0, None),
                3: ('O', 0, None), 4: ('O', 0, None)},
               {frozenset({0, 2}): (1, None), frozenset({1, 4}): (1, None),
                frozenset({2, 3}): (1, None), frozenset({2, 4}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.insertion(rct_cgrs, prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)


def test__reac__substitution():
    """ test graph.reac.substitution
    """
    # CH3OOH + CH3 => CH3OCH3 + OH
    rct_cgr = ({0: ('C', 3, None), 1: ('O', 1, None), 2: ('O', 0, None),
                3: ('C', 3, None)},
               {frozenset({0, 2}): (1, None), frozenset({1, 2}): (1, None)})
    prd_cgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('O', 0, None),
                3: ('O', 1, None)},
               {frozenset({0, 2}): (1, None), frozenset({1, 2}): (1, None)})

    rct_cgr = graph.explicit(rct_cgr)
    prd_cgr = graph.explicit(prd_cgr)

    rct_cgrs = graph.connected_components(rct_cgr)
    prd_cgrs = graph.connected_components(prd_cgr)

    tras, rct_idxs, prd_idxs = graph.reac.substitution(rct_cgrs, prd_cgrs)
    assert tras
    assert rct_idxs
    assert prd_idxs

    print("substitution")
    for tra in tras:
        print(tra)
        assert graph.backbone_isomorphic(
            graph.trans.apply(tra, rct_cgr), prd_cgr)


def test__prod__hydrogen_abstraction():
    """ test graph.reac.prod_hydrogen_abstraction
    """

    x_gra = automol.geom.graph(
        automol.inchi.geometry(
            automol.smiles.inchi('CC')))
    y_gra = automol.geom.graph(
        automol.inchi.geometry(
            automol.smiles.inchi('[OH]')))

    prod_gras = graph.reac.prod_hydrogen_abstraction(x_gra, y_gra)

    print('\n h abstraction')
    for gra in prod_gras:
        print(gra)


def test__prod__hydrogen_migration():
    """ test graph.reac.prod_hydrogen migration
    """

    gra = automol.geom.graph(
        automol.inchi.geometry(
            automol.smiles.inchi('CCCC[CH2]')))

    prod_gras = graph.reac.prod_hydrogen_migration(gra)

    print('\n h mig')
    for gra in prod_gras:
        print(gra)


def test__prod__addition():
    """ test graph.reac.prod_addition
    """

    x_gra = automol.geom.graph(
        automol.inchi.geometry(
            automol.smiles.inchi('C=C=C')))
    y_gra = automol.geom.graph(
        automol.inchi.geometry(
            automol.smiles.inchi('[H]')))

    prod_gras = graph.reac.prod_addition(x_gra, y_gra)

    print('\n addn')
    for gra in prod_gras:
        print(gra)


def test__prod__beta_scission():
    """ test graph.reac.prod_beta_scission
    """

    gra = automol.geom.graph(
        automol.inchi.geometry(
            automol.smiles.inchi('CCC[O]')))

    prod_gras = graph.reac.prod_beta_scission(gra)

    print('\n beta sci')
    for gra in prod_gras:
        print(gra)


def test__prod__homolytic_scission():
    """ test graph.reac.prod_homolytic_scission
    """

    gra = automol.geom.graph(
        automol.inchi.geometry(
            automol.smiles.inchi('CCCC')))

    prod_gras = graph.reac.prod_homolytic_scission(gra)

    print('\n homolyt scission')
    for gra in prod_gras:
        print(gra)


if __name__ == '__main__':
    # test__trans__is_stereo_compatible()
    # test__reac__hydrogen_migration()
    # test__reac__hydrogen_abstraction()
    # test__reac__addition()
    # test__reac__beta_scission()
    test__reac__elimination()
    # test__reac__insertion()
    # test__reac__substitution()
    # test__prod__hydrogen_abstraction()
    # test__prod__hydrogen_migration()
    # test__prod__addition()
    # test__prod__beta_scission()
    # test__prod__homolytic_scission()
