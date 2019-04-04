""" test the autoread.zmatrix module
"""
import autoread
import autoparse.pattern as app


def test_():
    """ test autoread.zmatrix
    """
    string = ('  Geometry (in Angstrom), charge = 0, multiplicity = 1:\n'
              '\n'
              '  O\n'
              '  O       1     R1\n'
              '  H       1     R2   2     A2\n'
              '  H       2     R2   1     A2   3     D3\n'
              '\n'
              '  A2    =  96.7725720000\n'
              '  D3    = 129.3669950000\n'
              '  R1    =  1.4470582953\n'
              '  R2    =  0.9760730000\n')

    syms, key_mat, name_mat, val_dct = autoread.zmatrix.read(
        string,
        start_ptt=(app.padded(app.escape('Geometry (in Angstrom),'),
                              app.NONNEWLINE)
                   + 2 * app.padded(app.NEWLINE)))
    assert syms == ('O', 'O', 'H', 'H')
    assert key_mat == ((None, None, None),
                       (1, None, None),
                       (1, 2, None),
                       (2, 1, 3))
    assert name_mat == ((None, None, None),
                        ('R1', None, None),
                        ('R2', 'A2', None),
                        ('R2', 'A2', 'D3'))
    assert val_dct == {
        'A2': 96.772572, 'D3': 129.366995, 'R1': 1.4470582953, 'R2': 0.976073}

    string = ('C \n'
              'O , 1 , R1 \n'
              'H , 1 , R2 , 2 , A2 \n'
              'H , 1 , R3 , 2 , A3 , 3 , D3 \n'
              'H , 1 , R4 , 2 , A4 , 3 , D4 \n'
              'H , 2 , R5 , 1 , A5 , 3 , D5 \n'
              '\n'
              'R1  = 2.67535  \n'
              'R2  = 2.06501   A2  = 109.528  \n'
              'R3  = 2.06501   A3  = 109.528   D3  = 120.808        \n'
              'R4  = 2.06458   A4  = 108.982   D4  = 240.404        \n'
              'R5  = 1.83748   A5  = 107.091   D5  = 299.596        \n')

    syms, key_mat, name_mat, val_dct = autoread.zmatrix.read(
        string,
        mat_entry_start_ptt=',',
        mat_entry_sep_ptt=',',
        setv_sep_ptt=app.padded(app.one_of_these(['', app.NEWLINE])))

    assert syms == ('C', 'O', 'H', 'H', 'H', 'H')
    assert key_mat == ((None, None, None),
                       (1, None, None),
                       (1, 2, None),
                       (1, 2, 3),
                       (1, 2, 3),
                       (2, 1, 3))
    assert name_mat == ((None, None, None),
                        ('R1', None, None),
                        ('R2', 'A2', None),
                        ('R3', 'A3', 'D3'),
                        ('R4', 'A4', 'D4'),
                        ('R5', 'A5', 'D5'))
    assert val_dct == {
        'R1': 2.67535, 'R2': 2.06501, 'A2': 109.528, 'R3': 2.06501,
        'A3': 109.528, 'D3': 120.808, 'R4': 2.06458, 'A4': 108.982,
        'D4': 240.404, 'R5': 1.83748, 'A5': 107.091, 'D5': 299.596}


def test__matrix():
    """ test autoread.zmatrix.matrix
    """
    string = (' comment:\n'
              ' --------\n'
              ' Symbolic Z-matrix:\n'
              ' Charge = 0 Multiplicity = 1\n'
              ' O\n'
              ' O     1 R1\n'
              ' H     1 R2  2 A2\n'
              ' H     2 R2  1 A2  3 D3  0\n'
              '  Variables:\n'
              ' R1     1.45405\n'
              '  Constants:\n'
              ' R2     0.97607\n'
              ' A2     96.77257\n'
              ' D3     129.367\n')
    syms, key_mat, name_mat = autoread.zmatrix.matrix.read(
        string,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('Symbolic Z-matrix:'), app.LINE, '']))

    assert syms == ('O', 'O', 'H', 'H')
    assert key_mat == ((None, None, None),
                       (1, None, None),
                       (1, 2, None),
                       (2, 1, 3))
    assert name_mat == ((None, None, None),
                        ('R1', None, None),
                        ('R2', 'A2', None),
                        ('R2', 'A2', 'D3'))

    string = (' -----------------------------------------------------------\n'
              '              Z-MATRIX (ANGSTROMS AND DEGREES)\n'
              '    nt   At 1      Lengt    2      Alph           Bet      J\n'
              ' -----------------------------------------------------------\n'
              '  1   1  C\n'
              '  2   2  O  1  1.454832( 1)\n'
              '  3   3  H  1  1.093067( 2)  2 111.219( 6)\n'
              '  4   4  H  1  1.090938( 3)  2 106.548( 7) 3 118.980( 10)  0\n'
              '  5   5  H  1  1.093052( 4)  2 111.241( 8) 3 238.017( 11)  0\n'
              '  6   6  O  2  1.360181( 5)  1 108.248( 9) 3  60.990( 12)  0\n'
              ' -----------------------------------------------------------')

    line_start_ptt = app.LINESPACES.join(2 * [app.UNSIGNED_INTEGER])
    entry_end_ptt = app.PADDING.join([
        app.escape('('), app.UNSIGNED_INTEGER, app.escape(')'),
        app.maybe(app.UNSIGNED_INTEGER)])

    syms, key_mat, name_mat = autoread.zmatrix.matrix.read(
        string,
        start_ptt=app.padded(app.NEWLINE).join([
            app.escape('Z-MATRIX (ANGSTROMS AND DEGREES)'),
            app.LINE, app.LINE, '']),
        name_ptt=app.FLOAT,
        line_start_ptt=line_start_ptt,
        entry_end_ptt=entry_end_ptt)

    assert syms == ('C', 'O', 'H', 'H', 'H', 'O')
    assert key_mat == ((None, None, None),
                       (1, None, None),
                       (1, 2, None),
                       (1, 2, 3),
                       (1, 2, 3),
                       (2, 1, 3))
    assert name_mat == ((None, None, None),
                        (1.454832, None, None),
                        (1.093067, 111.219, None),
                        (1.090938, 106.548, 118.98),
                        (1.093052, 111.241, 238.017),
                        (1.360181, 108.248, 60.99))
    print(name_mat)


def test__setval():
    """ test autoread.zmatrix.setval
    """
    string = ('    A2        =   96.7725720000\n'
              '    D3        =  129.3669950000\n'
              '    R1        =    1.4470582953\n'
              '    R2        =    0.9760730000\n')

    val_dct = autoread.zmatrix.setval.read(string)

    assert val_dct == {
        'A2': 96.772572, 'D3': 129.366995, 'R1': 1.4470582953, 'R2': 0.976073}

    string = ('              ----------------------------\n'
              '              !   Optimized Parameters   !\n'
              '              ! (Angstroms and Degrees)  !\n'
              ' -------------                            -----\n'
              ' !  Name     Value   Derivative information   !\n'
              ' ----------------------------------------------\n'
              ' !   R1     1.4057   -DE/DX =    0.0          !\n'
              ' !   R2     0.9761   -DE/DX =    0.0628       !\n'
              ' !   A2    96.7726   -DE/DX =    0.0552       !\n'
              ' !   D3   129.367    -DE/DX =    0.0019       !\n'
              ' ----------------------------------------------\n'
              ' GradGradGradGradGradGradGradGradGradGradGradGrad\n')

    start_ptt = app.padded(app.NEWLINE).join([
        app.escape('!   Optimized Parameters   !'),
        app.LINE, app.LINE, app.LINE, app.LINE, ''])

    val_dct = autoread.zmatrix.setval.read(
        string,
        start_ptt=start_ptt,
        entry_sep_ptt='',
        entry_start_ptt=app.escape('!'),
        sep_ptt=app.maybe(app.LINESPACES).join([
            app.escape('-DE/DX ='), app.FLOAT, app.escape('!'), app.NEWLINE]))

    assert val_dct == {
        'R1': 1.4057, 'R2': 0.9761, 'A2': 96.7726, 'D3': 129.367}

    string = ('R1  = 2.73454  \n'
              'R2  = 1.84451   A2  = 96.7726  \n'
              'R3  = 1.84451   A3  = 96.7726   D3  = 129.367        \n')

    val_dct = autoread.zmatrix.setval.read(
        string,
        sep_ptt=app.one_of_these(['', app.NEWLINE]))

    assert val_dct == {
        'R1': 2.73454, 'R2': 1.84451, 'A2': 96.7726, 'R3': 1.84451,
        'A3': 96.7726, 'D3': 129.367}


if __name__ == '__main__':
    # test__setval()
    # test_()
    test__matrix()
