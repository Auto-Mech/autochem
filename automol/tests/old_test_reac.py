""" test automol.reac
"""

import numpy
import automol


def test__stereo():
    """ test stereo functionality
    """

    # example 1
    rct_smis = ['FC=C(C(O)F)C(O)F', '[OH]']
    prd_smis = ['FC(O)[C](C(O)F)C(O)F']

    rxn_objs = automol.reac.with_structures_from_smiles(rct_smis, prd_smis)
    rxn, _, rct_geos, prd_geos = rxn_objs[0]

    # Complete stereo expansion for the reaction
    srxns = automol.reac.expand_stereo(rxn)
    print(len(srxns))
    print("Complete stereo expansion for the reaction:")
    srxn_chis = []
    for srxn in srxns:
        rct_gras = automol.reac.reactant_graphs(srxn)
        prd_gras = automol.reac.product_graphs(srxn)
        rct_chis = tuple(map(automol.graph.chi, rct_gras))
        prd_chis = tuple(map(automol.graph.chi, prd_gras))
        srxn_chis.append((rct_chis, prd_chis))
        print(rct_chis)
        print(prd_chis)
        print()

    assert set(srxn_chis) == {
        (('InChI=1S/C4H5F3O2/c5-1-2(3(6)8)4(7)9/h1,3-4,8-9H/'
          'b2-1-/t3-,4+/m0/s1', 'InChI=1S/HO/h1H'),
         ('InChI=1S/C4H6F3O3/c5-2(8)1(3(6)9)4(7)10/h2-4,8-10H/'
          't2-,3-,4+/m0/s1',)),
        (('InChI=1S/C4H5F3O2/c5-1-2(3(6)8)4(7)9/h1,3-4,8-9H/'
          'b2-1-/t3-,4+/m0/s1', 'InChI=1S/HO/h1H'),
         ('InChI=1S/C4H6F3O3/c5-2(8)1(3(6)9)4(7)10/h2-4,8-10H/'
          't2-,3-,4+/m1/s1',)),
        (('InChI=1S/C4H5F3O2/c5-1-2(3(6)8)4(7)9/h1,3-4,8-9H/'
          'b2-1-/t3-,4+/m1/s1', 'InChI=1S/HO/h1H'),
         ('InChI=1S/C4H6F3O3/c5-2(8)1(3(6)9)4(7)10/h2-4,8-10H/'
          't2-,3-,4+/m0/s1',)),
        (('InChI=1S/C4H5F3O2/c5-1-2(3(6)8)4(7)9/h1,3-4,8-9H/'
          'b2-1-/t3-,4+/m1/s1', 'InChI=1S/HO/h1H'),
         ('InChI=1S/C4H6F3O3/c5-2(8)1(3(6)9)4(7)10/h2-4,8-10H/'
          't2-,3-,4+/m1/s1',)),
        (('InChI=1S/C4H5F3O2/c5-1-2(3(6)8)4(7)9/h1,3-4,8-9H/'
          't3-,4-/m0/s1', 'InChI=1S/HO/h1H'),
         ('InChI=1S/C4H6F3O3/c5-2(8)1(3(6)9)4(7)10/h2-4,8-10H/'
          't2-,3-,4+/m0/s1',)),
        (('InChI=1S/C4H5F3O2/c5-1-2(3(6)8)4(7)9/h1,3-4,8-9H/'
          't3-,4-/m0/s1', 'InChI=1S/HO/h1H'),
         ('InChI=1S/C4H6F3O3/c5-2(8)1(3(6)9)4(7)10/h2-4,8-10H/'
          't2-,3-,4-/m0/s1',)),
        (('InChI=1S/C4H5F3O2/c5-1-2(3(6)8)4(7)9/h1,3-4,8-9H/'
          't3-,4-/m1/s1', 'InChI=1S/HO/h1H'),
         ('InChI=1S/C4H6F3O3/c5-2(8)1(3(6)9)4(7)10/h2-4,8-10H/'
          't2-,3-,4+/m1/s1',)),
        (('InChI=1S/C4H5F3O2/c5-1-2(3(6)8)4(7)9/h1,3-4,8-9H/'
          't3-,4-/m1/s1', 'InChI=1S/HO/h1H'),
         ('InChI=1S/C4H6F3O3/c5-2(8)1(3(6)9)4(7)10/h2-4,8-10H/'
          't2-,3-,4-/m1/s1',)),
    }

    # Assign reactant and product stereo from geometries.
    srxn = automol.reac.add_stereo_from_geometries(rxn, rct_geos, prd_geos)
    # Note that the original stereo assignments from the product geometries
    # could be inconsistent with the reactant stereo assignments.
    print('Consistent?', automol.reac.stereo_is_physical(srxn))

    # example 2
    rct_smis = ['FC=CC=CF', '[OH]']
    prd_smis = ['FC=C[CH]C(O)F']
    print("Reaction:", rct_smis, "=>", prd_smis)

    rxn_objs = automol.reac.with_structures_from_smiles(rct_smis, prd_smis)
    rxn, _, rct_geos, prd_geos = rxn_objs[0]

    # Complete stereo expansion for the reaction
    srxns = automol.reac.expand_stereo(rxn)
    print(len(srxns))
    print("Complete stereo expansion for the reaction:")
    print(rxn.forward_ts_graph)
    srxn_chis = []
    for srxn in srxns:
        rct_gras = automol.reac.reactant_graphs(srxn)
        prd_gras = automol.reac.product_graphs(srxn)
        rct_chis = tuple(map(automol.graph.chi, rct_gras))
        prd_chis = tuple(map(automol.graph.chi, prd_gras))
        srxn_chis.append((rct_chis, prd_chis))
        print(rct_chis)
        print(prd_chis)
        print("Physical?", automol.reac.stereo_is_physical(srxn))
        print()

    assert set(srxn_chis) == {
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1+,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1+,3-1+/t4-/m0/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1+,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1+,3-1+/t4-/m1/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1+,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1-,3-1+/t4-/m0/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1+,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1-,3-1+/t4-/m1/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1+,3-1+/t4-/m0/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1+,3-1+/t4-/m1/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1+,3-1-/t4-/m0/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1+,3-1-/t4-/m1/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1-,3-1+/t4-/m0/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1-,3-1+/t4-/m1/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1-,3-1-/t4-/m0/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2+', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1-,3-1-/t4-/m1/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2-', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1+,3-1-/t4-/m0/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2-', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1+,3-1-/t4-/m1/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2-', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1-,3-1-/t4-/m0/s1',)),
        (('InChI=1S/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1-,4-2-', 'InChI=1S/HO/h1H'),
         ('AMChI=1/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H/b2-1-,3-1-/t4-/m1/s1',)),
    }

    # Assign reactant and product stereo from geometries.
    srxn = automol.reac.add_stereo_from_geometries(rxn, rct_geos, prd_geos)
    # Note that the original stereo assignments from the product geometries
    # could be inconsistent with the reactant stereo assignments.
    print('Consistent?', automol.reac.stereo_is_physical(srxn))

    # example 3
    rct_smis = ['C(F)(Cl)-C(F)(Cl)O[O]']
    prd_smis = ['C(F)(Cl)=C(F)(Cl)', 'O[O]']

    rxn_objs = automol.reac.with_structures_from_smiles(rct_smis, prd_smis)
    rxn, _, rct_geos, prd_geos = rxn_objs[0]

    # Complete stereo expansion for the reaction
    srxns = automol.reac.expand_stereo(rxn)
    print(len(srxns))
    print("Complete stereo expansion for the reaction:")
    srxn_chis = []
    for srxn in srxns:
        rcts_gra = automol.reac.reactants_graph(srxn)
        prds_gra = automol.reac.products_graph(srxn)
        rcts_chi = automol.graph.chi(rcts_gra)
        prds_chi = automol.graph.chi(prds_gra)
        srxn_chis.append((rcts_chi, prds_chi))
        print(rcts_chi)
        print(prds_chi)
        print()

    assert set(srxn_chis) == {
        ('InChI=1S/C2HCl2F2O2/c3-1(5)2(4,6)8-7/h1H/t1-,2-/m0/s1',
         'InChI=1S/C2Cl2F2.HO2/c3-1(5)2(4)6;1-2/h;1H/b2-1-;'),
        ('InChI=1S/C2HCl2F2O2/c3-1(5)2(4,6)8-7/h1H/t1-,2+/m0/s1',
         'InChI=1S/C2Cl2F2.HO2/c3-1(5)2(4)6;1-2/h;1H/b2-1+;'),
        ('InChI=1S/C2HCl2F2O2/c3-1(5)2(4,6)8-7/h1H/t1-,2+/m1/s1',
         'InChI=1S/C2Cl2F2.HO2/c3-1(5)2(4)6;1-2/h;1H/b2-1+;'),
        ('InChI=1S/C2HCl2F2O2/c3-1(5)2(4,6)8-7/h1H/t1-,2-/m1/s1',
         'InChI=1S/C2Cl2F2.HO2/c3-1(5)2(4)6;1-2/h;1H/b2-1-;'),
    }

    # example 4
    rct_smis = ['N(F)-N(F)O[O]']
    prd_smis = ['N(F)=N(F)', 'O[O]']

    rxn_objs = automol.reac.with_structures_from_smiles(rct_smis, prd_smis)
    rxn, _, rct_geos, prd_geos = rxn_objs[0]

    # Complete stereo expansion for the reaction
    srxns = automol.reac.expand_stereo(rxn)
    print(len(srxns))
    print("Complete stereo expansion for the reaction:")
    srxn_chis = []
    for srxn in srxns:
        rcts_gra = automol.reac.reactants_graph(srxn)
        prds_gra = automol.reac.products_graph(srxn)
        rcts_chi = automol.graph.chi(rcts_gra)
        prds_chi = automol.graph.chi(prds_gra)
        srxn_chis.append((rcts_chi, prds_chi))
        print(rcts_chi)
        print(prds_chi)
        print()

    assert set(srxn_chis) == {
        ('AMChI=1/HF2N2O2/c1-3-4(2)6-5/h3H/t3-,4-/m0/s1',
         'InChI=1S/F2N2.HO2/c1-3-4-2;1-2/h;1H/b4-3+;'),
        ('AMChI=1/HF2N2O2/c1-3-4(2)6-5/h3H/t3-,4+/m0/s1',
         'InChI=1S/F2N2.HO2/c1-3-4-2;1-2/h;1H/b4-3-;'),
        ('AMChI=1/HF2N2O2/c1-3-4(2)6-5/h3H/t3-,4+/m1/s1',
         'InChI=1S/F2N2.HO2/c1-3-4-2;1-2/h;1H/b4-3-;'),
        ('AMChI=1/HF2N2O2/c1-3-4(2)6-5/h3H/t3-,4-/m1/s1',
         'InChI=1S/F2N2.HO2/c1-3-4-2;1-2/h;1H/b4-3+;'),
    }

    # example 5
    rct_smis = ['[CH2]CC=CC']
    prd_smis = ['[CH]=CC', 'C=C']

    rxn_objs = automol.reac.with_structures_from_smiles(rct_smis, prd_smis)
    rxn, _, rct_geos, prd_geos = rxn_objs[0]

    # Complete stereo expansion for the reaction
    srxns = automol.reac.expand_stereo(rxn)
    print(len(srxns))
    print("Complete stereo expansion for the reaction:")
    srxn_chis = []
    for srxn in srxns:
        rcts_gra = automol.reac.reactants_graph(srxn)
        prds_gra = automol.reac.products_graph(srxn)
        rcts_chi = automol.graph.chi(rcts_gra)
        prds_chi = automol.graph.chi(prds_gra)
        srxn_chis.append((rcts_chi, prds_chi))
        print(rcts_chi)
        print(prds_chi)
        print()

    assert set(srxn_chis) == {
        ('InChI=1S/C5H9/c1-3-5-4-2/h4-5H,1,3H2,2H3/b5-4-',
         'AMChI=1/C3H5.C2H4/c1-3-2;1-2/h1,3H,2H3;1-2H2/b3-1+;'),
        ('InChI=1S/C5H9/c1-3-5-4-2/h4-5H,1,3H2,2H3/b5-4+',
         'AMChI=1/C3H5.C2H4/c1-3-2;1-2/h1,3H,2H3;1-2H2/b3-1-;'),
    }


def test__add_stereo_from_unordered_geometries():
    """ test automol.reac.add_stereo_from_unordered_geometries()
    """
    rct_smis = ['C(F)(Cl)C(F)(Cl)', '[OH]']
    prd_smis = ['C(F)(Cl)[C](F)(Cl)', 'O']
    rxn = automol.reac.with_structures_from_smiles(rct_smis, prd_smis)[0][0]

    rct_geos = list(map(automol.inchi.geometry,
                        map(automol.smiles.inchi, rct_smis)))
    prd_geos = list(map(automol.inchi.geometry,
                        map(automol.smiles.inchi, prd_smis)))

    # mix the geometries up a bit so they no longer correspond
    def randomize_atom_ordering(geo):
        """ randomize atom ordering in a geometry
        """
        natms = automol.geom.count(geo)
        ord_dct = dict(enumerate(numpy.random.permutation(natms)))
        return automol.geom.reorder(geo, ord_dct)

    rct_geos = list(map(randomize_atom_ordering, reversed(rct_geos)))
    prd_geos = list(map(randomize_atom_ordering, reversed(prd_geos)))

    # This would break because the order doesn't match:
    # srxn = automol.reac.add_stereo_from_geometries(srxn, rct_geos, prd_geos)

    # We do this instead:
    srxn, order = automol.reac.add_stereo_from_unordered_geometries(
        rxn, rct_geos, prd_geos)
    if srxn is None:
        # The stereo must be inconsistent -- reflect coordinates for the
        # products (Only guaranteed to work for this particular reaction)
        prd_geos = list(map(automol.geom.reflect_coordinates, prd_geos))
        srxn, order = automol.reac.add_stereo_from_unordered_geometries(
            rxn, rct_geos, prd_geos)
    print(automol.reac.string(srxn))
    print(order)

    # Here's how we would reorder our geometries to match the reaction object:
    rct_order, prd_order = order
    rct_geos = [rct_geos[i] for i in rct_order]
    prd_geos = [prd_geos[i] for i in prd_order]

    # Now this should work:
    srxn = automol.reac.add_stereo_from_geometries(srxn, rct_geos, prd_geos)


def test__canonical_enantiomer():
    """ test reac.canonical_enantiomer
    """
    rct_smis = ['CC(OO)C(O[O])C(OO)C']
    prd_smis = ['CC(OO)C(OO)C(OO)[CH2]']

    rxn = automol.reac.with_structures_from_smiles(rct_smis, prd_smis)[0][0]

    # 2A. Full expansion -- includes non-canonical enantiomer reactions
    print("Full reaction expansion:")
    for srxn in automol.reac.expand_stereo(rxn, enant=True):
        rct_chis, prd_chis = automol.reac.chi(srxn)
        print(' +\n'.join(rct_chis) + " =>\n" + ' +\n'.join(prd_chis))

        # These functions operate directly on the reaction object:
        is_can = automol.reac.is_canonical_enantiomer(srxn)
        print(f"Canonical? {is_can}")
        # Convert it to a canonical enantiomer reaction like this
        srxn = automol.reac.canonical_enantiomer(srxn)
        assert automol.reac.is_canonical_enantiomer(srxn)

        # These are the equivalent functions for ChIs
        is_can = automol.chi.is_canonical_enantiomer_reaction(rct_chis,
                                                              prd_chis)
        print(f"Canonical? {is_can}")
        # Convert it to a canonical enantiomer reaction like this
        rct_chis, prd_chis = automol.chi.canonical_enantiomer_reaction(
            rct_chis, prd_chis)
        assert automol.chi.is_canonical_enantiomer_reaction(rct_chis, prd_chis)
        print()

    # 2B. Restricted expansion -- includes only canonical enantiomers
    print("Restricted reaction expansion:")
    for srxn in automol.reac.expand_stereo(rxn, enant=False):
        rct_chis, prd_chis = automol.reac.chi(srxn)
        print(' +\n'.join(rct_chis) + " =>\n" + ' +\n'.join(prd_chis))

        # Check canonicity for a reaction object
        assert automol.reac.is_canonical_enantiomer(srxn)

        # Check canonicity for reaction ChIs
        assert automol.chi.is_canonical_enantiomer_reaction(rct_chis, prd_chis)
        print()


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings("error")

    # test__add_stereo_from_unordered_geometries()
    # test__stereo()
    # test__canonical_enantiomer()
