""" test automol.automol.etrans
"""

import numpy
import automol


# Info for various target species
# Closed-Shell Hydrocarbons with different C-C p,s,t,q connectivities
CH1_INF = ('InChI=1S/C2H6/c1-2/h1-2H3', 1, 0)                     # CC
CH2_INF = ('InChI=1S/C8H18/c1-3-5-7-8-6-4-2/h3-8H2,1-2H3', 1, 0)  # CCCCCCCC
CYC_INF = ('InChI=1S/C6H12/c1-2-4-6-5-3-1/h1-6H2', 1, 0)          # C1CCCCC1
DMB_INF = ('InChI=1S/C6H14/c1-5-6(2,3)4/h5H2,1-4H3', 1, 0)        # CCC(C)(C)C
DMB2_INF = ('InChI=1S/C6H14/c1-5(2)6(3)4/h5-6H,1-4H3', 1, 0)      # CC(C)C(C)C
TRMB_INF = ('InChI=1S/C7H16/c1-6(2)7(3,4)5/h6H,1-5H3', 1, 0)
# ^ CC(C)C(C)(C)C
TEMB_INF = ('InChI=1S/C8H18/c1-7(2,3)8(4,5)6/h1-6H3', 1, 0)
# ^ CC(C)(C)C(C)(C)C
ISOPENT_INF = ('InChI=1S/C5H12/c1-4-5(2)3/h5H,4H2,1-3H3', 1, 0)  # CCC(C)C

# Open-Shell C,H and C,H,O species
RAD1_INF = ('InChI=1S/C3H5/c1-3-2/h3H2,1H3', 2, 0)               # CC[C]
RAD2_INF = ('InChI=1S/C3H7O2/c1-2-3-5-4/h2-3H2,1H3', 2, 0)       # CCCO[O]

# Species with various functional groups
COOH_INF = ('InChI=1S/C3H8O2/c1-2-3-5-4/h4H,2-3H2,1H3', 1, 0)    # CCCOO
ETH_INF = ('InChI=1S/C4H10O/c1-3-5-4-2/h3-4H2,1-2H3', 1, 0)      # CCOCC
COOH_ETH_INF = ('InChI=1S/C3H8O3/c1-5-2-3-6-4/h4H,2-3H2,1H3', 1, 0)  # COCCOO
ALC_INF = ('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 0, 1)                # CCO
CX_INF = ('InChI=1S/C4H9Cl/c1-2-3-4-5/h2-4H2,1H3', 1, 0)             # CCCCCl
EPOX_INF = ('InChI=1S/C7H14O/c1-3-4-5-7(2)6-8-7/h3-6H2,1-2H3', 1, 0)
# ^ CCCCC1(CO1)C

# Info for various baths (need to set a default of bath not viable
BATH_INF = ('InChI=1S/Ar', 0, 1)
BATH_INF2 = ('InChI=1S/N2/c1-2', 0, 1)


def test__model():
    """ test automol.etrans.effective_model
    """

    # Assess different targets for model selection
    model = automol.etrans.effective_model(CH1_INF[0], BATH_INF[0])
    assert model == (BATH_INF[0], 'n-alkane')

    model = automol.etrans.effective_model(RAD1_INF[0], BATH_INF[0])
    assert model == (BATH_INF[0], '1-alkyl')

    model = automol.etrans.effective_model(RAD2_INF[0], BATH_INF[0])
    assert model == (BATH_INF[0], '1-alkyl')

    model = automol.etrans.effective_model(COOH_INF[0], BATH_INF[0])
    assert model == (BATH_INF[0], 'n-hydroperoxide')

    model = automol.etrans.effective_model(ETH_INF[0], BATH_INF[0])
    assert model == (BATH_INF[0], 'ether')

    model = automol.etrans.effective_model(COOH_ETH_INF[0], BATH_INF[0])
    assert model == (BATH_INF[0], 'n-hydroperoxide')

    model = automol.etrans.effective_model(ALC_INF[0], BATH_INF[0])
    assert model == (BATH_INF[0], 'n-alcohol')

    model = automol.etrans.effective_model(EPOX_INF[0], BATH_INF[0])
    assert model == (BATH_INF[0], 'epoxide')

    # Halide groups have no model, just using alkanes
    model = automol.etrans.effective_model(CX_INF[0], BATH_INF[0])
    assert model == (BATH_INF[0], 'n-alkane')

    # Use a different bath gas
    model = automol.etrans.effective_model(CH1_INF[0], BATH_INF2[0])
    assert model == (BATH_INF2[0], 'n-alkane')


def test__effective_rotor():
    """ test effective rotor count
    """

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(CH1_INF[0]))
    assert numpy.isclose(2.0, n_eff)

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(CH2_INF[0]))
    assert numpy.isclose(8.0, n_eff)

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(TRMB_INF[0]))
    assert numpy.isclose(3.3333333, n_eff)

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(TEMB_INF[0]))
    assert numpy.isclose(3.0, n_eff)

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(DMB_INF[0]))
    assert numpy.isclose(3.3333333, n_eff)

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(DMB2_INF[0]))
    assert numpy.isclose(3.6666666, n_eff)

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(CYC_INF[0]))
    assert numpy.isclose(3.0, n_eff)

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(ISOPENT_INF[0]))
    assert numpy.isclose(4.0, n_eff)

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(COOH_INF[0]))
    assert numpy.isclose(3.666666, n_eff)

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(COOH_ETH_INF[0]))
    assert numpy.isclose(3.333333, n_eff)

    n_eff = automol.etrans.eff.effective_rotor_count(
        automol.inchi.geometry(ALC_INF[0]))
    assert numpy.isclose(2.333333, n_eff)


def test__estimate():
    """ test automol.etrans.eff
    """

    # Set values for calculation
    model = (BATH_INF[0], 'n-alcohol')
    n_eff = 2.333333333333

    tgt_geo = automol.inchi.geometry(ALC_INF[0])
    bath_geo = automol.inchi.geometry(BATH_INF[0])
    tgt_mass = automol.geom.total_mass(tgt_geo)
    bath_mass = automol.geom.total_mass(bath_geo)

    tgt_n_heavy = automol.geom.atom_count(tgt_geo, 'H', match=False)

    # Build the Lennard-Jones parameters
    sig, eps = automol.etrans.eff.lennard_jones_params(
        tgt_n_heavy, *model)
    ref_sig = 3.799479365827328
    ref_eps = 206.2796676062968
    assert numpy.isclose(sig, ref_sig)
    assert numpy.isclose(eps, ref_eps)

    # Build the alpha parameters
    edown_alpha, edown_n = automol.etrans.eff.alpha(
        n_eff, eps, sig, tgt_mass, bath_mass, *model)
    ref_edown_alpha = 138.20750436962405
    ref_edown_n = 0.6659477362988162
    assert numpy.isclose(edown_alpha, ref_edown_alpha)
    assert numpy.isclose(edown_n, ref_edown_n)


def test__estimate_no_model():
    """ test automol.etrans.eff for unsupported model
    """

    # Test model with bath gas where there are no models
    nheavy = 3
    model = ('InChI=1S/Kr', 'n-alkane')
    sig, eps = automol.etrans.eff.lennard_jones_params(nheavy, *model)
    assert sig is None and eps is None


def test__combine():
    """ test automol.etrans.combine.epsilon
        test automol.etrans.combine.sigma
    """

    ab_eps, bb_eps = 218.7, 156.3
    ab_sig, bb_sig = 3.6, 2.4

    ref_aa_eps = 306.01209213051817
    ref_aa_sig = 4.800000000000001
    aa_eps = automol.etrans.combine.epsilon(ab_eps, bb_eps)
    aa_sig = automol.etrans.combine.sigma(ab_sig, bb_sig)

    assert numpy.isclose(aa_eps, ref_aa_eps)
    assert numpy.isclose(aa_sig, ref_aa_sig)


test__effective_rotor()
