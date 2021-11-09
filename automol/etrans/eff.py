"""
Calculate an effective alpha parameter using certain rules
"""

import itertools
import numpy
from phydat import phycon
import automol.geom
import automol.inchi
import automol.graph
from automol.util import dict_
# from automol.etrans._par import LJ_DCT
from automol.etrans._par import LJ_EST_DCT
from automol.etrans._par import Z_ALPHA_EST_DCT
from automol.etrans._fxn import troe_lj_collision_frequency


# CALCULATE THE EFFECTIVE ALPHA VALUE
def alpha(n_eff, eps, sig, mass1, mass2, bath_model, tgt_model,
          empirical_factor=2.0):
    """ Calculate the alpha param using the method Jasper, et al.

        :param n_eff: number of effective rotors
        :type n_eff: int
        :p`aram zlj_dct: lennard-jones collision frequencies (cm-1?)
        :type zlj_dct: dict[float: float]
        :param bath_model: InChI string for bath gas species
        :type bath_model: str
        :param target_model: string denoting some class of target species
        :type target_model: str
        :param empirical_factor: correction for using 1DME versus 2DM2
        :type empirical_factor: float
    """

    # Calculate the Lennard-Jones frequencies
    red_mass = ((mass1 * mass2) / (mass1 + mass2)) * phycon.AMU2KG

    # Calculate Zalpha(Neff) at T = 300, 1000, 2000 K
    z_alphas_n_eff = _calculate_z_alpha_terms(n_eff, bath_model, tgt_model)

    # Calculate alpha = Zalpha(Neff) / Z(N) at T = 300, 1000, 2000 K
    # Empirical correction factor of (1/2) used for 1D Master Equations
    alpha_dct = {}
    for temp, z_alpha_n_eff in z_alphas_n_eff.items():
        zlj = troe_lj_collision_frequency(eps, sig, red_mass, temp)
        alpha_dct[temp] = (z_alpha_n_eff / zlj) / empirical_factor

    # Determine alpha and n for the e-down model
    edown_alpha, edown_n = _calculate_energy_down_exponent(alpha_dct)

    return edown_alpha, edown_n


def _calculate_z_alpha_terms(n_eff, bath_model, tgt_model):
    """ Calculate the [Z*alpha](N_eff)
    """

    def _z_alpha(coeff, n_eff):
        """ calculate an effective Z*alpha parameter
        """
        return ((coeff[0] * n_eff**(3) +
                 coeff[1] * n_eff**(2) +
                 coeff[2] * n_eff**(1) +
                 coeff[3]) / 1.0e9)

    # Need to put special values in for H2 here

    # Read the proper coefficients from the moldriver dct
    coeff_dct = dict_.values_in_multilevel_dct(
        Z_ALPHA_EST_DCT, bath_model, tgt_model)

    if coeff_dct is not None:
        # Calculate the three alpha terms
        z_alpha_dct = {}
        for temp, coeffs in coeff_dct.items():
            z_alpha_dct[temp] = _z_alpha(coeffs, n_eff)

    return z_alpha_dct


def _calculate_energy_down_exponent(alpha_dct):
    """ Calculate power n, for model: E_down = E_down_300 * (T/300)**n

        Does a least-squares for n to solve the linear equation
        ln(E_down/E_down_300) = [ln(T/300)] * n

        :param alpha_dct: temperature-dependent alpha parameters
        :type alpha_dct: dict[float: float]
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

    return edown_alpha, edown_n


# CALCULATE THE EFFECTIVE LENNARD-JONES SIGMA AND EPSILON
def lennard_jones_params(n_heavy, bath_model, tgt_model):
    """ Returns in angstrom and cm-1.

        :param n_heavy: Number of heavy atoms for a species
        :type n_heavy: int
        :param bath_model: InChI string for bath gas species
        :type bath_model: str
        :param target_model: string denoting some class of target species
        :type target_model: str
    """

    def _lj(param, n_heavy, expt):
        """ calculate an effective Lennard-Jones parameter
        """
        return param * n_heavy**(expt)

    # Need to put special values in for H2 here

    # Read the proper coefficients from the moldriver dct
    coeffs = dict_.values_in_multilevel_dct(
        LJ_EST_DCT, bath_model, tgt_model)

    if coeffs is not None:
        # Calculate the effective sigma and epsilon values
        sig = _lj(coeffs[0], n_heavy, coeffs[1])
        eps = _lj(coeffs[2], n_heavy, coeffs[3])
    else:
        sig, eps = None, None

    return sig, eps


# DETERMINE N_EFF USED FOR ALPHA AND LJ PARAM CALCULATIONS
def effective_rotor_count(geo):
    """ Calculate an effective N parameter using the given parametrization.

        :param geo: geometry (Bohr)
        :type geo: automol geometry data structure
        :rtype: float
    """

    # Conver the geo to a graph
    gra = automol.geom.graph(geo)
    symbs = automol.geom.symbols(geo)

    # Count the rotors
    (n_pp, n_ps, n_pt, n_pq,
     n_ss, n_st, n_sq,
     n_tt, n_tq,
     n_qq,
     n_co, n_oo,
     n_ss_ring, n_rings) = _rotor_counts(gra, symbs)

    # print('    - Rotor Counts for N_eff:')
    # print('       N_pp:{}, N_ps:{}, N_pt:{}, N_pq:{}'.format(
    #     n_pp, n_ps, n_pt, n_pq))
    # print('       N_ss:{}, N_st:{}, N_sq:{}'.format(n_ss, n_st, n_sq))
    # print('       N_tt:{}, N_tq:{}'.format(n_tt, n_tq))
    # print('       N_qq:{}'.format(n_qq))
    # print('       N_co:{}, N_oo:{}'.format(n_co, n_oo))
    # print('       N_ss_ring:{}, N_rings:{}'.format(n_ss_ring, n_rings))

    # Use the rotor counts and the coefficients to calculate Neff
    c_pp_ps_ss, c_pt_st, c_pq_sq = 1.0, 2.0/3.0, 1.0/3.0
    c_tt_tq_qq, c_co_oo, c_ss_ring = 0.0, 1.0/3.0, 1.0/2.0
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
    """ Count up various types of bonds for a structure.

        :param gra: molecular graph of species
        :type gra: automol graph data structure
        :param symbs: atomic symbols of species
        :type symbs: tuple(str)
        :rtype: tuple(float)
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
    ring_keys = set(itertools.chain(*automol.graph.rings_atom_keys(gra)))
    n_rings = len(rings)

    # Loop over the bonds and count the number of atoms
    neighbors = automol.graph.atoms_neighbor_atom_keys(gra)
    for bnd in automol.graph.bond_keys(gra):
        key1, key2 = bnd
        spair = (symbs[key1], symbs[key2])
        if spair == ('C', 'C'):
            # Figure out which neighbors are not hydrogen and count the number
            atom1_neighbors = neighbors[key1]
            numc1 = 0
            for neighbor1 in atom1_neighbors:
                if symbs[neighbor1] != 'H':
                    numc1 += 1
            atom2_neighbors = neighbors[key2]
            numc2 = 0
            for neighbor2 in atom2_neighbors:
                if symbs[neighbor2] != 'H':
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
    return (n_pp, n_ps, n_pt, n_pq,
            n_ss, n_st, n_sq,
            n_tt, n_tq,
            n_qq,
            n_co, n_oo,
            n_ss_ring, n_rings)
