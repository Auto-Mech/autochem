""" test automol.automol.etrans
"""

import numpy
import automol


# Info for various target species
CH_INF = ('InChI=1S/C8H18/c1-3-5-7-8-6-4-2/h3-8H2,1-2H3', 1, 0)  # CCCCCCCC
RAD1_INF = ('InChI=1S/C3H5/c1-3-2/h3H2,1H3', 2, 0)               # CC[C]
RAD2_INF = ('InChI=1S/C3H7O2/c1-2-3-5-4/h2-3H2,1H3', 2, 0)       # CCCO[O]
COOH_INF = ('InChI=1S/C3H8O2/c1-2-3-5-4/h4H,2-3H2,1H3', 1, 0)    # CCCOO
ETH_INF = ('InChI=1S/C4H10O/c1-3-5-4-2/h3-4H2,1-2H3', 1, 0)     # CCOCC
COOH_ETH_INF = ('InChI=1S/C3H8O3/c1-5-2-3-6-4/h4H,2-3H2,1H3', 1, 0)  # COCCOO
ALC_INF = ('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 0, 1)                # CCO
CX_INF = ('InChI=1S/C4H9Cl/c1-2-3-4-5/h2-4H2,1H3', 1, 0)             # CCCCCl
EPOX_INF = ('InChI=1S/C7H14O/c1-3-4-5-7(2)6-8-7/h3-6H2,1-2H3', 1, 0)
# ^ CCCCC1(CO1)C

# Info for various baths (need to set a default of bath not viable
BATH_INFO = ('InChI=1S/Ar', 0, 1)
BATH_INFO2 = ('InChI=1S/N2/c1-2', 0, 1)

TGT = 'n-alcohol'
BATH_INFO = ('InChI=1S/Ar', 0, 1)
TGT_GEO = automol.inchi.geometry(ALC_INF[0])
BATH_GEO = automol.inchi.geometry(BATH_INFO[0])
TGT_MASS = automol.geom.total_mass(TGT_GEO)
BATH_MASS = automol.geom.total_mass(BATH_GEO)
TGT_N_HEAVY = automol.geom.atom_count(TGT_GEO, 'H', match=False)
TEMPS = (300., 1000., 2000.)

# Fake LJ parameters to test combining
AB_EPS, BB_EPS = 218.7, 156.3
AB_SIG, BB_SIG = 3.6, 2.4


def __model():
    """ test automol.etrans.effective_model
    """

    # Assess different targets for model selection
    model = automol.etrans.effective_model(CH_INF, BATH_INFO)
    assert model == (BATH_INFO[0], 'n-alkane')

    model = automol.etrans.effective_model(RAD1_INF, BATH_INFO)
    assert model == (BATH_INFO[0], '1-alkyl')

    model = automol.etrans.effective_model(RAD2_INF, BATH_INFO)
    assert model == (BATH_INFO[0], '1-alkyl')

    model = automol.etrans.effective_model(COOH_INF, BATH_INFO)
    assert model == (BATH_INFO[0], 'n-hydroperoxide')

    model = automol.etrans.effective_model(ETH_INF, BATH_INFO)
    assert model == (BATH_INFO[0], 'ether')

    model = automol.etrans.effective_model(COOH_ETH_INF, BATH_INFO)
    assert model == (BATH_INFO[0], 'n-hydroperoxide')

    model = automol.etrans.effective_model(ALC_INF, BATH_INFO)
    assert model == (BATH_INFO[0], 'n-alcohol')
    n_eff = automol.etrans.eff.effective_rotor_count(TGT_GEO)
    ref_n_eff = 2.333333333333
    assert numpy.isclose(n_eff, ref_n_eff)

    model = automol.etrans.effective_model(EPOX_INF, BATH_INFO)
    assert model == (BATH_INFO[0], 'epoxide')

    model = automol.etrans.effective_model(CX_INF, BATH_INFO)
    assert model == (BATH_INFO[0], 'n-alkane')

    # Use a different bath gas
    model = automol.etrans.effective_model(CH_INF, BATH_INFO2)
    assert model == (BATH_INFO2[0], 'n-alkane')


def test__estimate():
    """ test automol.etrans.eff
    """

    # Set values for calculation
    model = (BATH_INFO[0], 'n-alcohol')
    n_eff = 2.333333333333

    # Build the Lennard-Jones parameters
    sig, eps = automol.etrans.eff.lennard_jones_params(
        TGT_N_HEAVY, *model)
    ref_sig = 3.799479365827328
    ref_eps = 206.2796676062968
    assert numpy.isclose(sig, ref_sig)
    assert numpy.isclose(eps, ref_eps)

    # Build the alpha parameters
    edown_alpha, edown_n = automol.etrans.eff.alpha(
        n_eff, eps, sig, TGT_MASS, BATH_MASS, *model)
    ref_edown_alpha = 138.20750436962405
    ref_edown_n = 0.6659477362988162
    assert numpy.isclose(edown_alpha, ref_edown_alpha)
    assert numpy.isclose(edown_n, ref_edown_n)


def test__combine():
    """ test automol.etrans.combine.epsilon
        test automol.etrans.combine.sigma
    """

    ref_aa_eps = 306.01209213051817
    ref_aa_sig = 4.800000000000001
    aa_eps = automol.etrans.combine.epsilon(AB_EPS, BB_EPS)
    aa_sig = automol.etrans.combine.sigma(AB_SIG, BB_SIG)

    assert numpy.isclose(aa_eps, ref_aa_eps)
    assert numpy.isclose(aa_sig, ref_aa_sig)
