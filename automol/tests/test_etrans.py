""" test automol.automol.etrans
"""

import numpy
import automol


# Geometries with different rotors
TGT_INFO = ('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 0, 1)
TGT = 'n-alcohol'
BATH_INFO = ('InChI=1S/Ar', 0, 1)
TGT_GEO = automol.inchi.geometry(TGT_INFO[0])
BATH_GEO = automol.inchi.geometry(BATH_INFO[0])
TGT_MASS = automol.geom.total_mass(TGT_GEO)
BATH_MASS = automol.geom.total_mass(BATH_GEO)
TGT_N_HEAVY = automol.geom.atom_count(TGT_GEO, 'H', match=False)
TEMPS = (300., 1000., 2000.)

# Fake LJ parameters to test combining
AB_EPS, BB_EPS = 218.7, 156.3
AB_SIG, BB_SIG = 3.6, 2.4


def test__estimate():
    """ test automol.etrans.eff
    """

    # Test the model determination
    model = automol.etrans.eff.determine_collision_model(TGT_INFO, BATH_INFO)
    assert model == (BATH_INFO[0], 'n-alcohol')

    # Build the effective N parameter
    n_eff = automol.etrans.eff.effective_rotor_count(TGT_GEO)
    ref_n_eff = 2.333333333333
    assert numpy.isclose(n_eff, ref_n_eff)

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
