""" pyx2z interface
"""
import pyx2z


def _x2z_molecule_from_geometry(geo):
    _mg = pyx2z.MolecGeom()
    for sym, xyz in geo:
        _atm = pyx2z.Atom(sym)
        _atm[0], _atm[1], _atm[2] = xyz
        _mg.push_back(_atm)
    _ps = pyx2z.PrimStruct(_mg)
    x2z_mol = pyx2z.MolecStruct(_ps)
    return x2z_mol


if __name__ == '__main__':
    GEO = (('C', (-0.70116587131, 0.0146227007587, -0.016166607003)),
           ('O', (1.7323365056, -0.9538524899, -0.5617192010)),
           ('H', (-0.9827048283, 0.061897979239, 2.02901783816)),
           ('H', (-0.8787925682, 1.91673409124, -0.80019507919)),
           ('H', (-2.12093033745, -1.21447973767, -0.87411360631)),
           ('H', (2.9512589894, 0.17507745634, 0.22317665541)))
    X2Z_MOL = _x2z_molecule_from_geometry(GEO)
    print(X2Z_MOL.zmatrix())
