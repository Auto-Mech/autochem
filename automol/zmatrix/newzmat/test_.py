""" test automol.zmatrix
"""

import automol
from automol.zmatrix.newzmat._bimol_ts import hydrogen_abstraction
from automol.zmatrix.newzmat._bimol_ts import addition
from automol.zmatrix.newzmat._bimol_ts import insertion
from automol.zmatrix.newzmat._bimol_ts import substitution
from automol.zmatrix.newzmat._unimol_ts import hydrogen_migration
from automol.zmatrix.newzmat._unimol_ts import beta_scission
from automol.zmatrix.newzmat._unimol_ts import concerted_unimol_elimination
from automol.zmatrix.newzmat._unimol_ts import ring_forming_scission
from automol.zmatrix.newzmat._util import shifted_standard_zmas_graphs


# ZMA Bank
C2H6_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('CC')))
C2H4_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('C=C')))
CH4_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('C')))
CH2_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('[CH2]')))
OH_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('[OH]')))
H_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('[H]')))
CH3_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('[CH3]')))
H2O_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('O')))
HO2_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('O[O]')))
CH2O_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('C=O')))
CH3CH2O_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('CC[O]')))
H2O2_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('OO')))
CH2COH_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('[CH2]CO')))
CH3CH2CH2CH2_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('CCC[CH2]')))
CH3CHCH2CH3_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('CC[CH]C')))
C3H8_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('CCC')))
CH3CH2OO_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('CCO[O]')))
CH2CH2OOH_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('[CH2]COO')))
cCH2OCH2_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('C1CO1')))


# BIMOL TS
def test__ts_hydrogen_abstraction():
    """ test zmatrix.ts.hydrogen_abstraction
    """

    rct_zmas = [CH4_ZMA, OH_ZMA]
    prd_zmas = [CH3_ZMA, H2O_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)

    tras, _, _, rtyp = automol.graph.reac.classify(rct_gras, prd_gras)
    print('\nrtyp', rtyp)

    zma_ret = hydrogen_abstraction(rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])


def test__ts_addition():
    """ test zmatrix.ts.addition
    """

    rct_zmas = [C2H4_ZMA, OH_ZMA]
    prd_zmas = [CH2COH_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)

    tras, _, _, rtyp = automol.graph.reac.classify(rct_gras, prd_gras)
    print('\nrtyp', rtyp)

    zma_ret = addition(rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])


def test__ts_insertion():
    """ test zmatrix.ts.insertion
    """

    rct_zmas = [C2H6_ZMA, CH2_ZMA]
    prd_zmas = [C3H8_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)

    tras, _, _, rtyp = automol.graph.reac.classify(rct_gras, prd_gras)
    print('\nrtyp', rtyp)

    zma_ret = insertion(rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])


def test__ts_substitution():
    """ test zmatrix.ts.substitution
    """

    rct_zmas = [H2O2_ZMA, H_ZMA]
    prd_zmas = [H2O_ZMA, OH_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)

    tras, _, _, rtyp = automol.graph.reac.classify(rct_gras, prd_gras)
    print('\nrtyp', rtyp)

    zma_ret = substitution(rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])


# UNIMOL TS
def test__ts_hydrogen_migration():
    """ test zmatrix.ts.hydrogen_migration
    """

    rct_zmas = [CH3CH2CH2CH2_ZMA]
    prd_zmas = [CH3CHCH2CH3_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    rct_zmas = [rct_zmas]
    prd_zmas = [prd_zmas]

    tras, _, _, rtyp = automol.graph.reac.classify(rct_gras, prd_gras)
    print('\nrtyp', rtyp)

    zma_ret = hydrogen_migration(rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])


def test__ts_beta_scission():
    """ test zmatrix.ts.beta_scission
    """

    rct_zmas = [CH3CH2O_ZMA]
    prd_zmas = [CH3_ZMA, CH2O_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)

    tras, _, _, rtyp = automol.graph.reac.classify(rct_gras, prd_gras)
    print('\nrtyp', rtyp)

    zma_ret = beta_scission(rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])


def test__ts_elimination():
    """ test zmatrix.ts.elimination
    """

    rct_zmas = [CH3CH2OO_ZMA]
    prd_zmas = [C2H4_ZMA, HO2_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    rct_zmas = [rct_zmas]
    prd_zmas = [prd_zmas]

    tras, _, _, rtyp = automol.graph.reac.classify(rct_gras, prd_gras)
    print('\nrtyp', rtyp)

    zma_ret = concerted_unimol_elimination(rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])


def test__ts_ring_forming_scission():
    """ test zmatrix.ts.ring_forming_scission
    """

    rct_zmas = [CH2CH2OOH_ZMA]
    prd_zmas = [cCH2OCH2_ZMA, OH_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    rct_zmas = [rct_zmas]
    prd_zmas = [prd_zmas]

    tras, _, _, rtyp = automol.graph.reac.classify(rct_gras, prd_gras)
    print('\nrtyp', rtyp)

    zma_ret = ring_forming_scission(rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])


if __name__ == '__main__':
    # BIMOL
    test__ts_hydrogen_abstraction()
    test__ts_addition()
    test__ts_substitution()
    # test__ts_insertion()
    # UNIMOL
    test__ts_hydrogen_migration()
    test__ts_beta_scission()
    test__ts_elimination()
    # test__ts_ring_forming_scission()
