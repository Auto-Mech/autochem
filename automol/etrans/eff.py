"""
Calculate an effective alpha parameter using certain rules
"""

import numpy
import automol
from automol.graph._util import ring_idxs
from automol.etrans._par import read_z_alpha_dct
from automol.etrans._par import read_lj_dct
from automol.etrans._fxn import troe_lj_collision_frequency


# CALCULATE THE EFFECTIVE ALPHA VALUE
def alpha(n_eff, zlj_dct, bath_model, tgt_model):
    """ Calculate the alpha param using the method in Ahren's paper
    """

    # Calculate Zalpha(Neff) at T = 300, 1000, 2000 K
    z_alphas_n_eff = _calc_z_alpha(n_eff, bath_model, tgt_model)

    # Calculate alpha = Zalpha(Neff) / Z(N) at T = 300, 1000, 2000 K
    # Empirical correction factor of (1/2) used for 1D Master Equations
    alphas = {}
    for temp, z_alpha_n_eff in z_alphas_n_eff.items():
        alphas[temp] = (z_alpha_n_eff / zlj_dct[temp]) / 2.0

    # Determine alpha and n for the e-down model
    edown_alpha, edown_n = _calc_edown_expt(alphas)

    return edown_alpha, edown_n


def _calc_z_alpha(n_eff, bath_model, tgt_model):
    """ Calculate the [Z*alpha](N_eff)
    """

    def _z_alpha(coeff, n_eff):
        """ calculate an effective Z*alpha parameter
        """
        return ((coeff[0] * n_eff**(3) +
                 coeff[1] * n_eff**(2) +
                 coeff[2] * n_eff**(1) +
                 coeff[3]) / 1.0e9)

    # Read the proper coefficients from the moldriver dct
    coeff_dct = read_z_alpha_dct(bath_model, tgt_model)

    # Calculate the three alpha terms
    z_alpha_dct = {}
    for temp, coeffs in coeff_dct.items():
        z_alpha_dct[temp] = _z_alpha(coeffs, n_eff)

    return z_alpha_dct


def lj_collision_frequency(sig, eps, mass1, mass2,
                           temps=(300., 1000., 2000.)):
    """ Calculate the collisin freq by Troe formula
    """

    # Calculate the reduced mass
    red_mass = ((mass1 * mass2) / (mass1 + mass2)) * 1.66053892173e-27
    # red_mass = ((mass1 * mass2) / (mass1 + mass2)) * phydat.phycon.AMU2KG

    #  Calculate the frequency at each temperature
    zlj_dct = {}
    for temp in temps:
        zlj_dct[temp] = troe_lj_collision_frequency(eps, sig, red_mass, temp)

    print('    - Collisional frequencies from LJ parameters')
    for temp, val in zlj_dct.items():
        print('       T = {0} K, zlj = {1:<.3e} cm3/s'.format(temp, val))

    return zlj_dct


def _calc_edown_expt(alpha_dct):
    """ Calculate power n, for model:
        E_down = E_down_300 * (T/300)**n

        Does a least-squares for n to solve the linear equation
        ln(E_down/E_down_300) = [ln(T/300)] * n
    """

    assert 300 in alpha_dct, (
        'Must have 300 K in alphas'
    )

    # Set the edown alpha to the value at 300 K
    edown_alpha = alpha_dct[300]

    # Build vectors and matrices used for the fitting
    temps = numpy.array(list(alpha_dct.keys()), dtype=numpy.float64)
    alphas = numpy.array(list(alpha_dct.values()), dtype=numpy.float64)

    n_vec = numpy.log(temps / 300.0)
    coeff_mat = numpy.array([n_vec], dtype=numpy.float64)
    coeff_mat = coeff_mat.transpose()

    edown_vec = numpy.log(alphas / edown_alpha)

    # Perform the least-squares fit
    theta = numpy.linalg.lstsq(coeff_mat, edown_vec, rcond=None)[0]

    # Set the the edown n value to the fitting parameter
    edown_n = theta[0]

    print('    - Alpha parameters from estimation')
    for temp, val in alpha_dct.items():
        print('       T = {0} K, alpha = {1:<.3f} cm-1'.format(temp, val))

    return edown_alpha, edown_n


# CALCULATE THE EFFECTIVE LENNARD-JONES SIGMA AND EPSILON
def lj_sig_eps(n_heavy, bath_model, tgt_model):
    """ Returns in angstrom and cm-1
    """

    def _lj(param, n_heavy, expt):
        """ calculate an effective Lennard-Jones parameter
        """
        return param * n_heavy**(expt)

    # Read the proper coefficients from the moldriver dct
    coeffs = read_lj_dct(bath_model, tgt_model)

    # Calculate the effective sigma and epsilon values
    sigma = _lj(coeffs[0], n_heavy, coeffs[1])
    epsilon = _lj(coeffs[2], n_heavy, coeffs[3])

    # Put in a check for zeros
    if sigma == 0 or epsilon == 0:
        sigma, epsilon = None, None

    return sigma, epsilon


# DETERMINE N_EFF USED FOR ALPHA AND LJ PARAM CALCULATIONS
def calc_n_eff(geo,
               c_pp_ps_ss=1.0, c_pt_st=(2.0/3.0), c_pq_sq=(1.0/3.0),
               c_tt_tq_qq=0.0, c_co_oo=(1.0/3.0), c_ss_ring=(1.0/2.0)):
    """ Calculate an effective N parameter using the given parametrization
    """

    # Conver the geo to a graph
    gra = automol.geom.graph(geo)
    symbs = automol.geom.symbols(geo)

    # Count the rotors
    [n_pp, n_ps, n_pt, n_pq,
     n_ss, n_st, n_sq,
     n_tt, n_tq,
     n_qq,
     n_co, n_oo,
     n_ss_ring, n_rings] = _rotor_counts(gra, symbs)

    print('    - Rotor Counts for N_eff:')
    print('       N_pp:{}, N_ps:{}, N_pt:{}, N_pq:{}'.format(
        n_pp, n_ps, n_pt, n_pq))
    print('       N_ss:{}, N_st:{}, N_sq:{}'.format(n_ss, n_st, n_sq))
    print('       N_tt:{}, N_tq:{}'.format(n_tt, n_tq))
    print('       N_qq:{}'.format(n_qq))
    print('       N_co:{}, N_oo:{}'.format(n_co, n_oo))
    print('       N_ss_ring:{}, N_rings:{}'.format(n_ss_ring, n_rings))

    # Use the rotor counts and the coefficients to calculate Neff
    n_eff = 1.0 + (
        c_pp_ps_ss * (n_pp + n_ps + n_ss) +
        c_pt_st * (n_pt + n_st) +
        c_pq_sq * (n_pq + n_sq) +
        c_tt_tq_qq * (n_tt + n_tq + n_qq) +
        c_co_oo * (n_co + n_oo) +
        c_ss_ring * n_ss_ring - n_rings
    )

    return n_eff


def _rotor_counts(gra, symbs):
    """ Count up various types of rotors
    """

    # Initialize the rotor counts
    n_pp, n_ps, n_pt, n_pq = 0, 0, 0, 0
    n_ss, n_st, n_sq = 0, 0, 0
    n_tt, n_tq = 0, 0
    n_qq = 0
    n_co, n_oo = 0, 0
    n_ss_ring, n_rings = 0, 0

    # Get the rings  and the number
    rings = automol.graph.rings(gra)
    ring_keys = set(ring_idxs(automol.graph.rings(gra)))
    n_rings = len(rings)

    # Loop over the bonds and count the number of atoms
    neighbors = automol.graph.atom_neighbor_keys(gra)
    for bnd in automol.graph.bond_keys(gra):
        key1, key2 = bnd
        spair = (symbs[key1], symbs[key2])
        if spair == ('C', 'C'):
            # Figure out which neighbors are carbons and count the number
            atom1_neighbors = neighbors[key1]
            numc1 = 0
            for neighbor1 in atom1_neighbors:
                if symbs[neighbor1] == 'C':
                    numc1 += 1
            atom2_neighbors = neighbors[key2]
            numc2 = 0
            for neighbor2 in atom2_neighbors:
                if symbs[neighbor2] == 'C':
                    numc2 += 1
            # Determine appropriate term to increment
            npair = (numc1, numc2)
            if npair == (1, 1):
                n_pp += 1
            elif npair in ((1, 2), (2, 1)):
                n_ps += 1
            elif npair in ((1, 3), (3, 1)):
                n_pt += 1
            elif npair in ((1, 4), (4, 1)):
                n_pq += 1
            elif npair == (2, 2):
                if {key1, key2} <= ring_keys:
                    n_ss_ring += 1
                else:
                    n_ss += 1
            elif npair in ((2, 3), (3, 2)):
                n_st += 1
            elif npair in ((2, 4), (4, 2)):
                n_sq += 1
            elif npair == (3, 3):
                n_tt += 1
            elif npair in ((3, 4), (4, 3)):
                n_tq += 1
            elif npair == (4, 4):
                n_qq += 1
        elif spair in (('C', 'O'), ('O', 'C')):
            n_co += 1
        elif spair == ('O', 'O'):
            n_oo += 1

    # Compile counts into a tuple
    return [n_pp, n_ps, n_pt, n_pq,
            n_ss, n_st, n_sq,
            n_tt, n_tq,
            n_qq,
            n_co, n_oo,
            n_ss_ring, n_rings]


# CHECKERS
BAD_ICHS = (
    'InChI=1S/H2/h1H'
)


def estimate_viable(well_ich, well_geo, bath_info):
    """ Assess whether we can estimate using the formula
    """

    # Initialize the models
    bath_model = bath_info[0]
    tgt_model = None

    # Build the graph
    well_gra = automol.geom.graph(well_geo)

    # Identify the the target model
    if well_ich not in BAD_ICHS:
        if automol.graph.radical_species(well_gra):
            tgt_model = '1-alkyl'
        elif automol.graph.hydrocarbon_species(well_gra):
            tgt_model = 'n-alkane'
        else:
            fgrp_dct = automol.graph.functional_group_dct(well_gra)
            if fgrp_dct[automol.graph.Fgroup.HYDROPEROXY]:
                tgt_model = 'n-hydroperoxide'
            elif fgrp_dct[automol.graph.Fgroup.EPOXIDE]:
                tgt_model = 'epoxide'
            elif fgrp_dct[automol.graph.Fgroup.ETHER]:
                tgt_model = 'ether'
            elif fgrp_dct[automol.graph.Fgroup.ALCOHOL]:
                tgt_model = 'n-alcohol'

        # For now, set model to alkanes if nothing found and set up return obj
        if tgt_model is None:
            tgt_model = 'n-alkane'

        ret = (bath_model, tgt_model)
    else:
        ret = None

    return ret
