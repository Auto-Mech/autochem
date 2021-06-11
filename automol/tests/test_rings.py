""" test ring functionality in graph
"""

from automol import graph
from automol import smiles
from automol import inchi
from automol import geom
from automol import zmat
import numpy

# cyclohexane
ICH1 = 'InChI=1S/C6H12/c1-2-4-6-5-3-1/h1-6H2'
# benzene
ICH2 = 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H'
# cyclic-ether
ICH3 = 'InChI=1S/C7H14O/c1-2-4-6-8-7-5-3-1/h1-7H2'
# 1-propylcyclopentane
ICH4 = 'InChI=1S/C8H16/c1-2-5-8-6-3-4-7-8/h8H,2-7H2,1H3'
# polycycle: cyclohexane+cyclopentane
ICH5 = 'InChI=1S/C9H16/c1-2-5-9-7-3-6-8(9)4-1/h8-9H,1-7H2/t8-,9-/m1/s1'

ZMA1 = geom.zmatrix(inchi.geometry(ICH1))
ZMA2 = geom.zmatrix(inchi.geometry(ICH2))
ZMA3 = geom.zmatrix(inchi.geometry(ICH3))
ZMA4 = geom.zmatrix(inchi.geometry(ICH4))
ZMA5 = geom.zmatrix(inchi.geometry(ICH5))


def test__rings():
    """ test graph.rings
    """
    c5h5n5o_cgr = (
        {0: ('C', 1, None), 1: ('C', 0, None), 2: ('C', 0, None),
         3: ('C', 0, None), 4: ('C', 0, None), 5: ('N', 2, None),
         6: ('N', 0, None), 7: ('N', 0, None), 8: ('N', 0, None),
         9: ('N', 1, None), 10: ('O', 1, None)},
        {frozenset({10, 4}): (1, None), frozenset({8, 2}): (1, None),
         frozenset({0, 6}): (1, None), frozenset({9, 3}): (1, None),
         frozenset({1, 2}): (1, None), frozenset({3, 7}): (1, None),
         frozenset({2, 5}): (1, None), frozenset({1, 6}): (1, None),
         frozenset({0, 7}): (1, None), frozenset({9, 4}): (1, None),
         frozenset({1, 3}): (1, None), frozenset({8, 4}): (1, None)})

    assert graph.rings(c5h5n5o_cgr) == (
        ({0: ('C', 1, None), 1: ('C', 0, None), 3: ('C', 0, None),
          6: ('N', 0, None), 7: ('N', 0, None)},
         {frozenset({0, 6}): (1, None), frozenset({3, 7}): (1, None),
          frozenset({0, 7}): (1, None), frozenset({1, 6}): (1, None),
          frozenset({1, 3}): (1, None)}),
        ({1: ('C', 0, None), 2: ('C', 0, None), 3: ('C', 0, None),
          4: ('C', 0, None), 8: ('N', 0, None), 9: ('N', 1, None)},
         {frozenset({8, 2}): (1, None), frozenset({9, 3}): (1, None),
          frozenset({1, 2}): (1, None), frozenset({9, 4}): (1, None),
          frozenset({1, 3}): (1, None), frozenset({8, 4}): (1, None)})
    )


def test__ring_systems():
    """ test graph.ring_systems
    """
    # molecule:
    # InChI=1S/C19H30/c1-2-4-14-10-12(9-13(14)3-1)5-7-17-16-8-6-15-11-
    # 18(16)19(15)17/h12-19H,1-11H2/
    gra = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 2, None),
            3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
            6: ('C', 2, None), 7: ('C', 2, None), 8: ('C', 1, None),
            9: ('C', 2, None), 10: ('C', 2, None), 11: ('C', 2, None),
            12: ('C', 1, None), 13: ('C', 1, None), 14: ('C', 2, None),
            15: ('C', 2, None), 16: ('C', 2, None), 17: ('C', 2, None),
            18: ('C', 2, None)},
           {frozenset({9, 13}): (1, None), frozenset({3, 6}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({11, 12}): (1, None),
            frozenset({13, 14}): (1, None), frozenset({3, 5}): (1, None),
            frozenset({0, 2}): (1, None), frozenset({1, 4}): (1, None),
            frozenset({12, 13}): (1, None), frozenset({0, 1}): (1, None),
            frozenset({1, 7}): (1, None), frozenset({12, 15}): (1, None),
            frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
            frozenset({16, 15}): (1, None), frozenset({4, 5}): (1, None),
            frozenset({16, 17}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({18, 4}): (1, None), frozenset({17, 14}): (1, None),
            frozenset({8, 10}): (1, None), frozenset({18, 10}): (1, None),
            frozenset({8, 11}): (1, None)})
    rsys = graph.ring_systems(gra)
    assert len(rsys) == 2

    rsy_rngs = list(map(graph.rings, rsys))
    assert tuple(map(len, rsy_rngs)) == (3, 2)


def test__ring_systems_decomposed_atom_keys():
    """ test graph.ring_systems_decomposed_atom_keys
    """
    # molecule:
    # InChI=1S/C19H30/c1-2-4-14-10-12(9-13(14)3-1)5-7-17-16-8-6-15-11-
    # 18(16)19(15)17/h12-19H,1-11H2/
    gra = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 2, None),
            3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
            6: ('C', 2, None), 7: ('C', 2, None), 8: ('C', 1, None),
            9: ('C', 2, None), 10: ('C', 2, None), 11: ('C', 2, None),
            12: ('C', 1, None), 13: ('C', 1, None), 14: ('C', 2, None),
            15: ('C', 2, None), 16: ('C', 2, None), 17: ('C', 2, None),
            18: ('C', 2, None)},
           {frozenset({9, 13}): (1, None), frozenset({3, 6}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({11, 12}): (1, None),
            frozenset({13, 14}): (1, None), frozenset({3, 5}): (1, None),
            frozenset({0, 2}): (1, None), frozenset({1, 4}): (1, None),
            frozenset({12, 13}): (1, None), frozenset({0, 1}): (1, None),
            frozenset({1, 7}): (1, None), frozenset({12, 15}): (1, None),
            frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
            frozenset({16, 15}): (1, None), frozenset({4, 5}): (1, None),
            frozenset({16, 17}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({18, 4}): (1, None), frozenset({17, 14}): (1, None),
            frozenset({8, 10}): (1, None), frozenset({18, 10}): (1, None),
            frozenset({8, 11}): (1, None)})

    decomps = graph.ring_systems_decomposed_atom_keys(gra)
    assert decomps == (((0, 1, 4, 5), (0, 2, 3, 5), (1, 7, 6, 3)),
                       ((8, 9, 13, 12, 11), (13, 14, 17, 16, 15, 12)))

# a1 = +/-q
# a2 = +/-a1


def test__ring_puckering():
    """ ring pucker
    """
    smi = 'CC1CCCCC1'
    ich = smiles.inchi(smi)
    geo = inchi.geometry(ich)
    zma = geom.zmatrix(geo)
    gra = zmat.graph(zma)
    rings_atoms = graph.rings_atom_keys(gra)
    val_dct = zmat.value_dictionary(zma)
    coos = zmat.coordinates(zma)
    geo = zmat.geometry(zma)
    da_names = zmat.dihedral_angle_names(zma)

    for ring_atoms in rings_atoms:
        rotate_hyds = []
        ngbs = graph.atom_sorted_neighbor_atom_keys(gra, ring_atoms[0])
        symbs = geom.symbols(geo)
        for ngb in ngbs:
            if symbs[ngb] == 'H':
                rotate_hyds.append(ngb)
        ring_value_dct = {}
        for da_name in da_names:
            da_idxs = list(coos[da_name])[0]
            if len(list(set(da_idxs) & set(ring_atoms))) == 4:
                print(da_name, da_idxs)
                ring_value_dct[da_name] = val_dct[da_name]
        dist_value_dct = {}
        for i, _ in enumerate(ring_atoms):
            dist_value_dct[i] = zmat.distance(
                zma, ring_atoms[i-1], ring_atoms[i])

        samp_range_dct = {}
        for key, value in ring_value_dct.items():
            samp_range_dct[key] = (value - numpy.pi/4, value + numpy.pi/4)

        print(zmat.samples(zma, 5, samp_range_dct))


def test__zmat_ring():
    """  test (add TS)
    """

    ref_rng_dct1 = {
        '1-2-5-8-11-14': {'D7': [0.16168073524433701, 1.7324770620392336],
                          'D10': [4.550708773750596, 6.121505100545493],
                          'D13': [0.16167968490856266, 1.7324760117034592]}
    }
    ref_rng_dct2 = {
        '1-2-4-6-8-10': {'D5': [5.497778966967656, 7.068575293762552],
                         'D7': [-0.7853916433944105, 0.7854046834004861],
                         'D9': [-0.7853965979951805, 0.7853997287997161]}
    }
    ref_rng_dct3 = {
        '1-2-5-8-11-14-15-18': {'D7': [4.107068942834604, 5.677865269629501],
                                'D10': [0.2823332978102957, 1.853129624605192],
                                'D13': [3.7747759383709774, 5.345572265165874],
                                'D14': [1.3922420091755172, 2.963038335970413],
                                'D17': [4.0726811920433565, 5.643477518838253]}
    }
    ref_rng_dct4 = {
        '1-2-5-8-11': {'D7': [4.965369609998152, 6.536165936793049],
                       'D10': [-0.06067570281684087, 1.5101206239780556]}
    }
    ref_rng_dct5 = {
        '1-2-5-8-11': {'D7': [4.877800020154778, 6.4485963469496745],
                       'D10': [-0.012877983213259836, 1.5579183435816368]},
        '5-8-21-18-15-9': {'D17': [-0.33350277464453076, 1.2372935521503658],
                           'D20': [4.488064445840583, 6.058860772635479]}}

    # Get lists of atoms in the ring
    rng_atoms1 = zmat.all_rings_atoms(ZMA1, zrxn=None)
    rng_atoms2 = zmat.all_rings_atoms(ZMA2, zrxn=None)
    rng_atoms3 = zmat.all_rings_atoms(ZMA3, zrxn=None)
    rng_atoms4 = zmat.all_rings_atoms(ZMA4, zrxn=None)
    rng_atoms5 = zmat.all_rings_atoms(ZMA5, zrxn=None)

    # Sampling ranges (includes dihedral calls)
    rng_dct1 = zmat.all_rings_dct(ZMA1, rng_atoms1)
    rng_dct2 = zmat.all_rings_dct(ZMA2, rng_atoms2)
    rng_dct3 = zmat.all_rings_dct(ZMA3, rng_atoms3)
    rng_dct4 = zmat.all_rings_dct(ZMA4, rng_atoms4)
    rng_dct5 = zmat.all_rings_dct(ZMA5, rng_atoms5)

    # Check distances (includes distance calc)
    # still need a dist check failure for testing
    rng_chk1 = zmat.all_rings_distances_reasonable(ZMA1, rng_atoms1)
    rng_chk2 = zmat.all_rings_distances_reasonable(ZMA2, rng_atoms2)
    rng_chk3 = zmat.all_rings_distances_reasonable(ZMA3, rng_atoms3)
    rng_chk4 = zmat.all_rings_distances_reasonable(ZMA4, rng_atoms4)
    rng_chk5 = zmat.all_rings_distances_reasonable(ZMA5, rng_atoms5)


# TODO
# def __fragment_ring():
#     """ test
#     """
#     # automol.geom.fragment_ring_geo(geo)
#     pass
# def __ring_check():
#     """ test
#     """
#     # automol.geom.ring_angles_passes(geo, ring_atoms, thresh=ATHRESH):


if __name__ == '__main__':
    # test__rings()
    # test__ring_systems()
    # test__ring_systems_decomposed_atom_keys()
    # test__ring_puckering()
    test__zmat_ring()
