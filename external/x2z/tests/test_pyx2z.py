""" test the pyx2z module
"""
import numpy
import pyx2z


def test__AtomBase_mass():
    """ test pyx2z.AtomBase.mass()
    """
    assert numpy.allclose(pyx2z.AtomBase('C').mass(), 21874.661832)
    assert numpy.allclose(pyx2z.AtomBase('C', 13).mass(), 23703.6661089)


def test__Atom___getitem__():
    """ test pyx2z.Atom.__getitem__() and __setitem()
    """
    a = pyx2z.Atom('C', 13)
    a[0] = 1.0
    a[1] = 1.1
    a[2] = 1.2
    assert numpy.allclose(a[0], 1.0)
    assert numpy.allclose(a[1], 1.1)
    assert numpy.allclose(a[2], 1.2)


def test__MolecGeom_size():
    """ test pyx2z.MolecGeom.size()
    """
    m = pyx2z.MolecGeom()
    assert m.size() == 0


def test__MolecGeom_push_back():
    """ test pyx2z.MolecGeom.push_back()
    """
    m = pyx2z.MolecGeom()
    m.push_back(pyx2z.Atom('C'))
    assert m.size() == 1


def test__MolecOrient_sym_num():
    """ test pyx2z.MolecOrient.sym_num()
    """
    asymbs = ['O', 'H', 'H']
    coords = [(-1.2516025626,  2.3683550357,  0.0000000000),
              (-0.2816025626,  2.3683550357,  0.0000000000),
              (-1.5749323743,  3.2380324089, -0.2828764736)]
    m = _molec_geom_obj(asymbs, coords)
    o = pyx2z.MolecOrient(m)
    assert o.sym_num() == 2


def test__MolecOrient_is_enantiomer():
    """ test pyx2z.MolecOrient.is_enantiomer()
    """
    asymbs = ['O', 'H', 'H']
    coords = [(-1.2516025626, 2.3683550357, 0.0000000000),
              (-0.2816025626, 2.3683550357, 0.0000000000),
              (-1.5749323743, 3.2380324089, -0.2828764736)]
    m = _molec_geom_obj(asymbs, coords)
    o = pyx2z.MolecOrient(m)
    assert o.is_enantiomer() is False


def test__MolecOrient_is_plane():
    """ test pyx2z.MolecOrient.is_plane()
    """
    asymbs = ['O', 'H', 'H']
    coords = [(-1.2516025626, 2.3683550357, 0.0000000000),
              (-0.2816025626, 2.3683550357, 0.0000000000),
              (-1.5749323743, 3.2380324089, -0.2828764736)]
    m = _molec_geom_obj(asymbs, coords)
    o = pyx2z.MolecOrient(m)
    assert o.is_plane() is True


def test__MolecOrient_is_linear():
    """ test pyx2z.MolecOrient.is_linear()
    """
    asymbs = ['O', 'H', 'H']
    coords = [(-1.2516025626, 2.3683550357, 0.0000000000),
              (-0.2816025626, 2.3683550357, 0.0000000000),
              (-1.5749323743, 3.2380324089, -0.2828764736)]
    m = _molec_geom_obj(asymbs, coords)
    o = pyx2z.MolecOrient(m)
    assert o.is_linear() is False


def test__MolecOrient_size():
    """ test pyx2z.MolecOrient.size()
    """
    asymbs = ['O', 'H', 'H']
    coords = [(-1.2516025626, 2.3683550357, 0.0000000000),
              (-0.2816025626, 2.3683550357, 0.0000000000),
              (-1.5749323743, 3.2380324089, -0.2828764736)]
    m = _molec_geom_obj(asymbs, coords)
    o = pyx2z.MolecOrient(m)
    assert o.size() == 3


def test__PrimStruct_is_connected():
    """ test pyx2z.PrimStruct.is_connected()
    """
    asymbs = ['O', 'H', 'H']
    coords = [(-1.2516025626, 2.3683550357, 0.0000000000),
              (-0.2816025626, 2.3683550357, 0.0000000000),
              (-1.5749323743, 3.2380324089, -0.2828764736)]
    m = _molec_geom_obj(asymbs, coords)
    p = pyx2z.PrimStruct(m)
    assert p.is_connected(0, 0) is True
    assert p.is_connected(0, 1) is True
    assert p.is_connected(0, 2) is True
    assert p.is_connected(1, 0) is True
    assert p.is_connected(1, 1) is True
    assert p.is_connected(1, 2) is False
    assert p.is_connected(2, 0) is True
    assert p.is_connected(2, 1) is False
    assert p.is_connected(2, 2) is True


def test__PrimStruct_connected_group():
    """ test pyx2z.PrimStruct.connected_group()
    """
    asymbs = ['O', 'H', 'H']
    coords = [(-1.2516025626, 2.3683550357, 0.0000000000),
              (-0.2816025626, 2.3683550357, 0.0000000000),
              (-1.5749323743, 3.2380324089, -0.2828764736)]
    m = _molec_geom_obj(asymbs, coords)
    p = pyx2z.PrimStruct(m)
    assert p.connected_group() == [[2, 0, 1]]


def test__PrimStruct_group_stoicheometry():
    """ test pyx2z.PrimStruct.group_stoicheometry()
    """
    asymbs = ['C', 'O', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H']
    coords = [(-2.2704178657, 2.9586871732, -0.0000000000),
              (-2.0230332171, 1.7537576097, 0.0000000000),
              (-1.3903791691, 3.9451234997, 0.7188849720),
              (-0.5605839137, 3.4182459749, 1.1980823377),
              (-1.9710320856, 4.4623547554, 1.4866301464),
              (-0.9857357677, 4.6646475866, 0.0028769813),
              (-3.4679340657, 3.5185797384, -0.7188849720),
              (-4.0230068872, 2.7073743052, -1.1980823377),
              (-3.1380484986, 4.2227540733, -1.4866301464),
              (-4.1233452839, 4.0204635187, -0.0028769813)]
    m = _molec_geom_obj(asymbs, coords)
    p = pyx2z.PrimStruct(m)
    assert p.group_stoicheometry([5, 2, 3, 4]) == 'C1H3'
    assert p.group_stoicheometry([9, 6, 7, 8]) == 'C1H3'
    assert p.group_stoicheometry([5, 2, 0, 1, 3, 4]) == 'C2H3O1'
    assert p.group_stoicheometry([9, 6, 0, 1, 7, 8]) == 'C2H3O1'


def test__MolecStruct_size():
    """ test pyx2z.MolecStruct.size()
    """
    asymbs = ['O', 'H', 'H']
    coords = [(-1.2516025626, 2.3683550357, 0.0000000000),
              (-0.2816025626, 2.3683550357, 0.0000000000),
              (-1.5749323743, 3.2380324089, -0.2828764736)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    assert s.size() == 3
    asymbs = ['H']
    coords = [(0., 0., 0.)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    assert s.size() == 1


def test__MolecStruct_resonance_averaged_bond_order():
    """ test pyx2z.MolecStruct.resonance_averaged_bond_order()
    """
    asymbs = ['O', 'H', 'H']
    coords = [(-1.2516025626, 2.3683550357, 0.0000000000),
              (-0.2816025626, 2.3683550357, 0.0000000000),
              (-1.5749323743, 3.2380324089, -0.2828764736)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    assert s.resonance_averaged_bond_order(0, 1) == 1.
    assert s.resonance_averaged_bond_order(0, 2) == 1.
    assert s.resonance_averaged_bond_order(1, 2) == 0.


def test__MolecStruct_is_radical():
    """ test pyx2z.MolecStruct.is_radical()
    """
    asymbs = ['O', 'H']
    coords = [(-1.2516025626, 2.3683550357, 0.0000000000),
              (-1.5749323743, 3.2380324089, -0.2828764736)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    assert s.is_radical(0) is True
    assert s.is_radical(1) is False
    asymbs = ['H']
    coords = [(0., 0., 0.)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    assert s.is_radical(0) is True


def test__MolecStruct_rotation_bond():
    """ test pyx2z.MolecStruct.rotation_bond()
    """
    asymbs = ['C', 'O', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H']
    coords = [(-2.2704178657, 2.9586871732, -0.0000000000),
              (-2.0230332171, 1.7537576097, 0.0000000000),
              (-1.3903791691, 3.9451234997, 0.7188849720),
              (-0.5605839137, 3.4182459749, 1.1980823377),
              (-1.9710320856, 4.4623547554, 1.4866301464),
              (-0.9857357677, 4.6646475866, 0.0028769813),
              (-3.4679340657, 3.5185797384, -0.7188849720),
              (-4.0230068872, 2.7073743052, -1.1980823377),
              (-3.1380484986, 4.2227540733, -1.4866301464),
              (-4.1233452839, 4.0204635187, -0.0028769813)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    assert s.rotation_bond() == {4: [[9, 6, 0, 1, 7, 8], [5, 2, 3, 4]],
                                 7: [[9, 6, 7, 8], [5, 2, 0, 1, 3, 4]]}


def test__MolecStruct_atom_ordering():
    """ test pyx2z.MolecStruct.atom_ordering()
    """
    asymbs = ['C', 'O', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H']
    coords = [(-2.2704178657, 2.9586871732, -0.0000000000),
              (-2.0230332171, 1.7537576097, 0.0000000000),
              (-1.3903791691, 3.9451234997, 0.7188849720),
              (-0.5605839137, 3.4182459749, 1.1980823377),
              (-1.9710320856, 4.4623547554, 1.4866301464),
              (-0.9857357677, 4.6646475866, 0.0028769813),
              (-3.4679340657, 3.5185797384, -0.7188849720),
              (-4.0230068872, 2.7073743052, -1.1980823377),
              (-3.1380484986, 4.2227540733, -1.4866301464),
              (-4.1233452839, 4.0204635187, -0.0028769813)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    assert s.atom_ordering() == [0, 1, 2, 6, 3, 4, 5, 7, 8, 9]


def test__MolecStruct_resonance_count():
    """ test pyx2z.MolecStruct.resonance_count()
    """
    asymbs = ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'H']
    coords = [(1.10206, 0.05263, 0.02517),
              (2.44012, 0.03045, 0.01354),
              (3.23570, 0.06292, 1.20436),
              (2.86296, -0.38925, 2.11637),
              (4.29058, 0.30031, 1.12619),
              (0.54568, -0.01805, -0.90370),
              (0.53167, 0.14904, 0.94292),
              (2.97493, -0.03212, -0.93001)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    assert s.resonance_count() == 2


def test__MolecStruct_bond_order():
    """ test pyx2z.MolecStruct.bond_order()
    """
    asymbs = ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'H']
    coords = [(1.10206, 0.05263, 0.02517),
              (2.44012, 0.03045, 0.01354),
              (3.23570, 0.06292, 1.20436),
              (2.86296, -0.38925, 2.11637),
              (4.29058, 0.30031, 1.12619),
              (0.54568, -0.01805, -0.90370),
              (0.53167, 0.14904, 0.94292),
              (2.97493, -0.03212, -0.93001)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    assert s.bond_order(0, 1, 0) == 2
    assert s.bond_order(1, 2, 0) == 1
    assert s.bond_order(0, 1, 1) == 1
    assert s.bond_order(1, 2, 1) == 2


def test__zmatrix_string():
    """ test pyx2z.zmatrix_sring()
    """
    asymbs = ['C', 'O', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H']
    coords = [(-2.2704178657, 2.9586871732, -0.0000000000),
              (-2.0230332171, 1.7537576097, 0.0000000000),
              (-1.3903791691, 3.9451234997, 0.7188849720),
              (-0.5605839137, 3.4182459749, 1.1980823377),
              (-1.9710320856, 4.4623547554, 1.4866301464),
              (-0.9857357677, 4.6646475866, 0.0028769813),
              (-3.4679340657, 3.5185797384, -0.7188849720),
              (-4.0230068872, 2.7073743052, -1.1980823377),
              (-3.1380484986, 4.2227540733, -1.4866301464),
              (-4.1233452839, 4.0204635187, -0.0028769813)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    string = pyx2z.zmatrix_string(s)
    print(string)


def test__rotational_bond_coordinates():
    """ test pyx2z.rotational_bond_coordinates()
    """
    asymbs = ['C', 'O', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H']
    coords = [(-2.2704178657, 2.9586871732, -0.0000000000),
              (-2.0230332171, 1.7537576097, 0.0000000000),
              (-1.3903791691, 3.9451234997, 0.7188849720),
              (-0.5605839137, 3.4182459749, 1.1980823377),
              (-1.9710320856, 4.4623547554, 1.4866301464),
              (-0.9857357677, 4.6646475866, 0.0028769813),
              (-3.4679340657, 3.5185797384, -0.7188849720),
              (-4.0230068872, 2.7073743052, -1.1980823377),
              (-3.1380484986, 4.2227540733, -1.4866301464),
              (-4.1233452839, 4.0204635187, -0.0028769813)]
    m = _molec_geom_obj(asymbs, coords)
    s = pyx2z.MolecStruct(pyx2z.PrimStruct(m))
    coords = pyx2z.rotational_bond_coordinates(s)
    assert coords == ['D4', 'D7']


def _molec_geom_obj(asymbs, coords):
    _mg = pyx2z.MolecGeom()
    for asymb, xyz in zip(asymbs, coords):
        _a = _atom_obj(asymb, xyz)
        _mg.push_back(_a)
    return _mg


def _atom_obj(asymb, xyz):
    _ang2bohr = 1.8897259886
    _a = pyx2z.Atom(asymb)
    _a[0], _a[1], _a[2] = numpy.multiply(xyz, _ang2bohr)
    return _a


if __name__ == '__main__':
    test__MolecStruct_resonance_count()
    test__MolecStruct_bond_order()
    test__MolecStruct_size()
    test__MolecStruct_is_radical()
    test__zmatrix_string()
    test__rotational_bond_coordinates()
