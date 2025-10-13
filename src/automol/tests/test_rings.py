""" test ring functionality in graph
tests on new functions in geom and zmat
"""

import os
from pathlib import Path
import numpy
from automol import geom, graph, inchi, smiles, zmat

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

GEO1 = inchi.geometry(ICH1)
GEO2 = inchi.geometry(ICH2)
GEO3 = inchi.geometry(ICH3)
GEO4 = inchi.geometry(ICH4)
GEO5 = inchi.geometry(ICH5)
ZMA1 = geom.zmatrix(GEO1)
ZMA2 = geom.zmatrix(GEO2)
ZMA3 = geom.zmatrix(GEO3)
ZMA4 = geom.zmatrix(GEO4)
ZMA5 = geom.zmatrix(GEO5)

PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, "data")


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


def __zmat_ring():
    """  test (add TS)
    """

    def _chk_ring_dct(ring_dct, ref_ring_dct):
        """ Ring dictionaries by checking the keys and subkeys and tha
            the floats match in the arrays.
        """
        for key, rkey in zip(ring_dct.keys(), ref_ring_dct.keys()):
            assert key == rkey
            rdct, ref_dct = ring_dct[key], ref_ring_dct[rkey]
            for key2, rkey2 in zip(rdct.keys(), ref_dct.keys()):
                assert key2 == rkey2
                assert numpy.allclose(rdct[key2], ref_dct[rkey2],
                                      atol=0.0001, rtol=0.0)

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
                           'D20': [4.488064445840583, 6.058860772635479]}
    }

    # Get lists of atoms in the ring
    rng_atoms1 = zmat.all_rings_atoms(ZMA1, tsg=None)
    rng_atoms2 = zmat.all_rings_atoms(ZMA2, tsg=None)
    rng_atoms3 = zmat.all_rings_atoms(ZMA3, tsg=None)
    rng_atoms4 = zmat.all_rings_atoms(ZMA4, tsg=None)
    rng_atoms5 = zmat.all_rings_atoms(ZMA5, tsg=None)

    # Sampling ranges (includes dihedral calls)
    rng_dct1 = zmat.all_rings_dct(ZMA1, rng_atoms1)
    rng_dct2 = zmat.all_rings_dct(ZMA2, rng_atoms2)
    rng_dct3 = zmat.all_rings_dct(ZMA3, rng_atoms3)
    rng_dct4 = zmat.all_rings_dct(ZMA4, rng_atoms4)
    rng_dct5 = zmat.all_rings_dct(ZMA5, rng_atoms5)

    _chk_ring_dct(rng_dct1, ref_rng_dct1)
    _chk_ring_dct(rng_dct2, ref_rng_dct2)
    _chk_ring_dct(rng_dct3, ref_rng_dct3)
    _chk_ring_dct(rng_dct4, ref_rng_dct4)
    _chk_ring_dct(rng_dct5, ref_rng_dct5)

    # Check distances (includes distance calc)
    # still need a dist check failure for testing
    assert zmat.all_rings_distances_reasonable(ZMA1, rng_atoms1)
    assert zmat.all_rings_distances_reasonable(ZMA2, rng_atoms2)
    assert zmat.all_rings_distances_reasonable(ZMA3, rng_atoms3)
    assert zmat.all_rings_distances_reasonable(ZMA4, rng_atoms4)
    assert zmat.all_rings_distances_reasonable(ZMA5, rng_atoms5)


def __geom_ring():
    """ test
    """

    # Check fragments
    ref_frag1 = (
        ('C', (-4.678309211005585, 0.9507582225493368, 1.3283943630808774)),
        ('C', (-5.121751782500965, -0.938130104907219, -0.8105253645907926)),
        ('C', (-1.8183184542355333, 1.124802798603796, 1.679473859905406)),
        ('H', (-5.596045421215773, 0.3417726406642164, 3.080387789063853)),
        ('H', (-5.472315138734762, 2.795412023268741, 0.8240154865962166)),
        ('C', (-2.5431899261233295, -2.0739380012437003, -1.3930195681295054)),
        ('H', (-5.852499087712589, 0.05531532655919256, -2.475545119599078)),
        ('H', (-6.503039672160678, -2.388190468828352, -0.2948252443802334)),
        ('C', (-0.7061241675458045, 0.061496234243248286, -0.764480412567453)),
        ('H', (-2.394903577788051, -2.688224956361526, -3.3620755359797467)),
        ('H', (-2.2076515854285743, -3.721079537831081, -0.17992394099843007)),
        ('C', (2.0244139428425885, -0.8586191572573748, -0.5328743916778614)),
        ('H', (-0.8215175624787503, 1.5073535945698067, -2.249882522830806)),
        ('H', (-1.2484554395000416, -0.024690004182670054, 3.307711745513383)),
        ('H', (-1.213466703672248, 3.0695743548570653, 2.0397346958319655)))
    ref_frag2 = (
        ('C', (-4.75275159254756, -0.8859368224809074, 0.05204272094591431)),
        ('C', (-2.6316145805140074, -2.781457031057775, 0.617675595973341)),
        ('C', (-3.469373724017983, 1.644203318184678, -0.5508463526763928)),
        ('H', (-6.020017911796573, -0.6793342485754604, 1.6757584497553977)),
        ('H', (-5.894535909176483, -1.541597262442682, -1.545311764072632)),
        ('C', (-0.2860088109004335, -1.148940006016929, 0.961220215463222)),
        ('H', (-3.0413520136248446, -3.922314862417004, 2.293255681709448)),
        ('H', (-2.4141209511227446, -4.0740160313602765, -0.9873049704822309)),
        ('C', (2.2490731922647558, -2.457979419683331, 0.5439361499900985)),
        ('C', (-0.7009866699077159, 0.9649583926305777, -0.9504147216562486)),
        ('H', (-0.3133961766730456, -0.3543813002328926, 2.882320427181015)),
        ('C', (1.200527146395269, 3.101000420082161, -0.5817249641001333)),
        ('H', (-0.4809624718722801, 0.2063149270378129, -2.873628003146178)),
        ('H', (-4.297962624106833, 2.5542573283204937, -2.212847874618906)),
        ('H', (-3.6981964276399215, 2.9346249752146134, 1.054281230413613)),
        ('C', (4.277907565282426, -0.6204973771653151, -0.42814819707015134)),
        ('C', (3.775171063722607, 2.1216311210353376, 0.3399582461004284)),
        ('H', (4.349010775297676, -0.7282474510540757, -2.4981877728699335)),
        ('H', (6.139118315609488, -1.2222378958530145, 0.25232279110138667)),
        ('H', (3.8464710144147634, 2.255925498026601, 2.408439361912304)),
        ('H', (5.279790586767728, 3.3507823389200784, -0.3773631616284865)),
        ('H', (1.4381199118411645, 4.1005566214301385, -2.3820194923439555)),
        ('H', (0.5015487139955048, 4.508608348000816, 0.7687730045421981)),
        ('H', (2.075693517927493, -4.028187328133698, -0.7970212520506802)),
        ('H', (2.8688480603817905, -3.297736252410089, 2.3348346516273333)))

    frag1 = geom.ring_fragments_geometry(GEO4)
    frag2 = geom.ring_fragments_geometry(GEO5)

    assert geom.almost_equal_dist_matrix(frag1, ref_frag1)
    assert geom.almost_equal_dist_matrix(frag2, ref_frag2)

    # Check angles passing
    # rng_atoms1 = zmat.all_rings_atoms(ZMA1, zrxn=None)
    # rng_atoms2 = zmat.all_rings_atoms(ZMA2, zrxn=None)
    # rng_atoms3 = zmat.all_rings_atoms(ZMA3, zrxn=None)
    # rng_atoms4 = zmat.all_rings_atoms(ZMA4, zrxn=None)
    # rng_atoms5 = zmat.all_rings_atoms(ZMA5, zrxn=None)
    # assert geom.all_rings_angles_reasonable(GEO1, rng_atoms1)
    # assert geom.all_rings_angles_reasonable(GEO2, rng_atoms2)
    # assert geom.all_rings_angles_reasonable(GEO3, rng_atoms3)
    # assert geom.all_rings_angles_reasonable(GEO4, rng_atoms4)
    # assert geom.all_rings_angles_reasonable(GEO5, rng_atoms5)
    # # Make fake ring with a bad angle
    # bad_geo = (
    #     ('C',  (-5.5877937580, -0.9968691886, -0.4989724332)),
    #     ('C',  (-2.7396917422, -0.5869879179, -0.3495489041)),
    #     ('H',  (-6.2131071344, -2.2170574292, 1.0990444734)),
    #     ('H',  (-6.2375188836, -1.8516712887, -2.3026059381)),
    #     ('C',  (-5.5653523999, 1.1009499352, 0.6068939227)),
    #     ('H',  (-1.6185907638, -2.3610282320, -0.3017615543)),
    #     ('C',  (-2.5439758909, 0.9534327614, 2.0676109788)),
    #     ('H',  (-2.1053962916, 0.5765661495, -1.9857880223)),
    #     ('C',  (-4.8046407570, 2.7722969017, 1.9610434891)),
    #     ('H',  (-7.3189290824, 1.5474187146, 1.4538721624)),
    #     ('H',  (-5.6843643376, 2.5317510457, -1.5414484200)),
    #     ('H',  (-0.6915593935, 1.8574930579, 2.4592560449)),
    #     ('H',  (-2.9167136073, -0.4107171297, 3.6299854071)),
    #     ('H',  (-5.7060177389, 2.8094515351, 3.8587924532)),
    #     ('H',  (-4.2677394999, 4.7261055569, 1.4049815516)))
    # rng_atoms = ((0, 1, 4, 6, 8),)
    # assert not geom.all_rings_angles_reasonable(bad_geo, rng_atoms)


def test__cremerpople():
    ring_geos = [GEO1,GEO2,GEO3,GEO4]
    for geo in ring_geos:
        ring_subgeo = geom.ring_only_geometry(geo)
        coord = [cord for _,cord in ring_subgeo]
        crem_pop,z = geom.cremer_pople_params(coord)
        q,phi = crem_pop
        assert len(z) == len(ring_subgeo)
        assert len(q)+len(phi) == len(ring_subgeo)-3


def test__dbscan_clustering():
    # chair_geo = (
    #             ('C',   (-0.498449,   1.373130,  -0.231044)),
    #             ('C',   (-1.438792,   0.254970,   0.230618)),
    #             ('H',   (-0.518909,   1.429388,  -1.328912)),
    #             ('H',   (-0.850553,   2.342729,   0.140418)),
    #             ('C',   (-0.940111,  -1.118642,  -0.230303)),
    #             ('H',   (-2.454323,   0.434869,  -0.141376)),
    #             ('H',   (-1.497800,   0.265878,   1.328485)),
    #             ('C',   ( 0.498874,  -1.373373,   0.230436)),
    #             ('H',   (-1.603521,  -1.908152,   0.142029)),
    #             ('H',   (-0.979109,  -1.164870,  -1.328153)),
    #             ('C',   ( 1.438724,  -0.254593,  -0.230841)),
    #             ('H',   ( 0.519731,  -1.430001,   1.328298)),
    #             ('H',   ( 0.850748,  -2.342753,  -0.141769)),
    #             ('C',   ( 0.939869,   1.118473,   0.231205)),
    #             ('H',   ( 2.454439,  -0.434609,   0.140627)),
    #             ('H',   ( 1.497365,  -0.264640,  -1.328740)),
    #             ('H',   ( 1.603264,   1.908325,  -0.140437)),
    #             ('H',   ( 0.977975,   1.164049,   1.329104)))
    # boat_geo = (
    #             ('C',   (-1.522738,   0.002041,  -0.001183)),
    #             ('C',   (-0.661994,  -1.219595,  -0.385613)),
    #             ('H',   (-2.177579,   0.263471,  -0.840807)),
    #             ('H',   (-2.181223,  -0.256410,   0.836540)),
    #             ('C',   ( 0.658915,  -1.220294,   0.386819)),
    #             ('H',   (-0.440951,  -1.197783,  -1.460997)),
    #             ('H',   (-1.218479,  -2.145790,  -0.205260)),
    #             ('C',   ( 1.522632,  -0.002036,  -0.001022)),
    #             ('H',   ( 0.438393,  -1.195663,   1.462232)),
    #             ('H',   ( 1.213086,  -2.148226,   0.208202)),
    #             ('C',   ( 0.661923,   1.219541,  -0.385839)),
    #             ('H',   ( 2.181557,   0.256728,   0.836251)),
    #             ('H',   ( 2.177168,  -0.263487,  -0.840920)),
    #             ('C',   (-0.658856,   1.220286,   0.386858)),
    #             ('H',   ( 1.218438,   2.145690,  -0.205311)),
    #             ('H',   ( 0.441008,   1.197793,  -1.461228)),
    #             ('H',   (-1.212746,   2.148550,   0.208977)),
    #             ('H',   (-0.437966,   1.195471,   1.462199)))
    with open(os.path.join(DAT_PATH, "c6h12_ring_samples.xyz"), "r") as f:
        traj = f.read()
    geos = [geoi for geoi,_ in geom.from_xyz_trajectory_string(traj)]

    rings_atoms = graph.rings_atom_keys(geom.graph(geos[0]))
    unique_geos = geom.dbscan(geos,rings_atoms,2.0,min_samples=1)

    assert len(unique_geos) == 5

    # I have not added check to checks_with_crest as it would require a CREST run


if __name__ == "__main__":
   # test__ring_puckering()
    test__cremerpople()
   # test__dbscan_clustering()
