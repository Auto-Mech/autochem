""" test ring functionality in graph
"""
from automol import graph
from automol import smiles
from automol import inchi
from automol import geom
from automol import zmat
import numpy


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


if __name__ == '__main__':
    # test__rings()
    # test__ring_systems()
    # test__ring_systems_decomposed_atom_keys()
    test__ring_puckering()
