""" test the autoread module
"""
import numpy
import autoread
import autoparse.pattern as app


def test__energy():
    """ test autoread.energy
    """
    string = (' @RKS Final Energy:  -149.44268093948725\n'
              '\n'
              ' => Energetics <=\n'
              '\n'
              ' Nuclear Repulsion Energy =   36.4973050010413047\n'
              ' One-Electron Energy =      -279.4237697165179384\n'
              ' Two-Electron Energy =       108.0155180600336564\n'
              ' Total Energy =             -149.4426809394872464\n')

    ene = autoread.energy.read(string, app.escape('@RKS Final Energy:'))
    assert ene == -149.44268093948725

    string = (' No special actions if energy rises.\n'
              ' SCF Done: E(RB3LYP) =  -149.442473330  A.U. after 9 cycles\n'
              '           NFock=  9 Conv=0.14D-08     -V/T= 2.0124\n'
              ' QCSCF skips out because SCF is already converged.\n')
    ene = autoread.energy.read(string, app.escape('E(RB3LYP) ='))
    assert ene == -149.442473330


def test__geom():
    """ test autoread.geom
    """
    string = ('                  Standard orientation:\n'
              ' --------------------------------------------------------\n'
              ' Center  Atomic  Atomic         Coordinates (Angstroms)\n'
              ' Number  Number   Type         X          Y          Z\n'
              ' --------------------------------------------------------\n'
              '      1       8       0   -0.000000   0.723527  -0.046053\n'
              '      2       8       0    0.000000  -0.723527  -0.046053\n'
              '      3       1       0    0.876174   0.838634   0.368421\n'
              '      4       1       0   -0.876174  -0.838634   0.368421\n'
              ' --------------------------------------------------------\n')

    nums, xyzs = autoread.geom.read(
        string,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('Standard orientation:'),
            app.LINE, app.LINE, app.LINE, app.LINE, '']),
        sym_ptt=app.UNSIGNED_INTEGER,
        line_start_ptt=app.UNSIGNED_INTEGER,
        line_sep_ptt=app.UNSIGNED_INTEGER,)

    assert nums == (8, 8, 1, 1)
    assert xyzs == ((-0.0, 0.723527, -0.046053),
                    (0.0, -0.723527, -0.046053),
                    (0.876174, 0.838634, 0.368421),
                    (-0.876174, -0.838634, 0.368421))

    string = (' Final (previous) structure:\n'
              ' Cartesian Geometry (in Angstrom)\n'
              '     O     0.0000000000   0.7028389815   0.0245676525\n'
              '     O     0.0000000000  -0.7028389815   0.0245676525\n'
              '     H    -0.8761735478   0.8179459165  -0.3899064740\n'
              '     H     0.8761735478  -0.8179459165  -0.3899064740\n')

    syms, xyzs = autoread.geom.read(
        string,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('Final (previous) structure:'), app.LINE, '']))

    assert syms == ('O', 'O', 'H', 'H')
    assert xyzs == ((0.0, 0.7028389815, 0.0245676525),
                    (0.0, -0.7028389815, 0.0245676525),
                    (-0.8761735478, 0.8179459165, -0.389906474),
                    (0.8761735478, -0.8179459165, -0.389906474))


def test__matrix():
    """ test autoread.matrix
    """
    # gaussian gradient
    string = (' ***** Axes restored to original set *****\n'
              ' ------------------------------------------------------------\n'
              ' Center   Atomic              Forces (Hartrees/Bohr)\n'
              ' Number   Number         X              Y              Z\n'
              ' ------------------------------------------------------------\n'
              '    1       8       0.000000000   -0.000000000   -0.061240635\n'
              '    2       1       0.000000000   -0.021691602    0.030622057\n'
              '    3       1      -0.000000000    0.021691602    0.030622057\n'
              ' ------------------------------------------------------------\n'
              ' Cartesian Forces:  Max     0.061240635 RMS     0.027012111\n')

    mat = autoread.matrix.read(
        string,
        start_ptt=app.padded(app.NEWLINE).join([
            app.padded(app.escape('Forces (Hartrees/Bohr)'), app.NONNEWLINE),
            app.LINE, app.LINE, '']),
        line_start_ptt=app.LINESPACES.join([app.UNSIGNED_INTEGER] * 2))

    assert numpy.allclose(mat, ((0.0, -0.0, -0.061240635),
                                (0.0, -0.021691602, 0.030622057),
                                (-0.0, 0.021691602, 0.030622057)))

    # gaussian hessian
    string = (' The second derivative matrix:\n'
              '             X1    Y1    Z1    X2    Y2\n'
              '      X1   -0.21406\n'
              '      Y1   -0.00000  2.05336\n'
              '      Z1   -0.00000  0.12105  0.19177\n'
              '      X2   -0.06169 -0.00000  0.00000  0.03160\n'
              '      Y2    0.00000 -0.09598 -0.05579 -0.00000  0.12501\n'
              '      Z2    0.00000  0.08316 -0.38831 -0.00000 -0.06487\n'
              '      X3    0.27574  0.00000  0.00000  0.03009 -0.00000\n'
              '      Y3    0.00000 -1.95737 -0.06525  0.00000 -0.02902\n'
              '      Z3    0.00000 -0.20421  0.19654 -0.00000  0.12066\n'
              '             Z2    X3    Y3    Z3\n'
              '      Z2    0.44623\n'
              '      X3    0.00000 -0.30583\n'
              '      Y3   -0.01829 -0.00000  1.98640\n'
              '      Z3   -0.05792 -0.00000  0.08354 -0.13862\n'
              ' ITU= 0\n'
              '   Eigenvalues ---  0.06664  0.66895  3.25121\n')

    comp_ptt = app.one_of_these(['X', 'Y', 'Z']) + app.UNSIGNED_INTEGER
    mat = autoread.matrix.read(
        string,
        start_ptt=(app.escape('The second derivative matrix:') +
                   app.lpadded(app.NEWLINE)),
        block_start_ptt=(app.series(comp_ptt, app.LINESPACES) +
                         app.padded(app.NEWLINE)),
        line_start_ptt=comp_ptt,
        tril=True)

    assert numpy.allclose(
        mat,
        ((-0.21406, 0., 0., -0.06169, 0., 0., 0.27574, 0., 0.),
         (0., 2.05336, 0.12105, 0., -0.09598, 0.08316, 0., -1.95737, -0.20421),
         (0., 0.12105, 0.19177, 0., -0.05579, -0.38831, 0., -0.06525, 0.19654),
         (-0.06169, 0., 0., 0.0316, 0., 0., 0.03009, 0., 0.),
         (0., -0.09598, -0.05579, 0., 0.12501, -0.06487, 0., -0.02902,
          0.12066),
         (0., 0.08316, -0.38831, 0., -0.06487, 0.44623, 0., -0.01829,
          -0.05792),
         (0.27574, 0., 0., 0.03009, 0., 0., -0.30583, 0., 0.),
         (0., -1.95737, -0.06525, 0., -0.02902, -0.01829, 0., 1.9864,
          0.08354),
         (0., -0.20421, 0.19654, 0., 0.12066, -0.05792, 0., 0.08354,
          -0.13862)))

    # psi4 gradient
    string = ('  ## Gradient (Symmetry 0) ##\n'
              '  Irrep: 1 Size: 3 x 3\n'
              '\n'
              '            1         2         3\n'
              '\n'
              '    1     0.000     0.000     0.997\n'
              '    2     0.000    -0.749    -0.488\n'
              '    3    -0.000     0.749    -0.488\n')

    mat = autoread.matrix.read(
        string,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('## Gradient (Symmetry 0) ##'),
            app.LINE, '', app.LINE, '', '']),
        line_start_ptt=app.UNSIGNED_INTEGER)

    assert numpy.allclose(
        mat, ((0.0, 0.0, 0.997), (0.0, -0.749, -0.488), (-0.0, 0.749, -0.488)))

    # psi4 hessian
    string = ('-------------------------------------------\n'
              ' ## Hessian (Symmetry 0) ##\n'
              ' Irrep: 1 Size: 9 x 9\n'
              '\n'
              '       1     2     3     4     5 \n'
              '\n'
              '  1   0.000   0.000   0.000   0.000   0.000\n'
              '  2   0.000   0.959   0.000   0.000  -0.452\n'
              '  3   0.000   0.000   0.371   0.000   0.222\n'
              '  4   0.000   0.000   0.000   0.000   0.000\n'
              '  5   0.000  -0.479   0.279   0.000   0.455\n'
              '  6   0.000   0.251  -0.185   0.000  -0.247\n'
              '  7   0.000   0.000   0.000   0.000   0.000\n'
              '  8   0.000  -0.479  -0.279   0.000  -0.003\n'
              '  9   0.000  -0.251  -0.185   0.000   0.025\n'
              '\n'
              '       6     7     8     9\n'
              '\n'
              '  1   0.000   0.000   0.000   0.000\n'
              '  2   0.519   0.000  -0.477  -0.230\n'
              '  3  -0.555   0.000  -0.279  -0.128\n'
              '  4   0.000   0.000   0.000   0.000\n'
              '  5  -0.256   0.000  -0.017   0.051\n'
              '  6   0.607   0.000  -0.012   0.090\n'
              '  7   0.000   0.000   0.000   0.000\n'
              '  8  -0.263   0.000   0.494   0.279\n'
              '  9   0.947   0.000   0.292   0.137\n')

    mat = autoread.matrix.read(
        string,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('## Hessian (Symmetry 0) ##'), app.LINE, '']),
        block_start_ptt=app.padded(app.NEWLINE).join([
            '', app.series(app.UNSIGNED_INTEGER, app.LINESPACES), '', '']),
        line_start_ptt=app.UNSIGNED_INTEGER)

    assert numpy.allclose(
        mat,
        ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
         (0.0, 0.959, 0.0, 0.0, -0.452, 0.519, 0.0, -0.477, -0.23),
         (0.0, 0.0, 0.371, 0.0, 0.222, -0.555, 0.0, -0.279, -0.128),
         (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
         (0.0, -0.479, 0.279, 0.0, 0.455, -0.256, 0.0, -0.017, 0.051),
         (0.0, 0.251, -0.185, 0.0, -0.247, 0.607, 0.0, -0.012, 0.09),
         (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
         (0.0, -0.479, -0.279, 0.0, -0.003, -0.263, 0.0, 0.494, 0.279),
         (0.0, -0.251, -0.185, 0.0, 0.025, 0.947, 0.0, 0.292, 0.137)))


if __name__ == '__main__':
    # test__energy()
    # test__geom()
    test__matrix()
