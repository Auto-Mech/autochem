""" Use various estimation formula of Jasper to determine
    energy transfer parameters: sigma, epsilon, and alpha
"""

import itertools

import numpy

from phydat import phycon

from .. import chi, geom, graph
from ..graph import FunctionalGroup
from ._fxn import troe_lj_collision_frequency
from ._par import (
    D0_GRP_LST,
    LJ_DCT,
    LJ_EST_DCT,
    Z_ALPHA_EST_DCT,
    ZROT_DCT,
)


# CALCULATE THE EFFECTIVE ALPHA VALUE
def alpha(n_eff, eps, sig, mass1, mass2, collider_set, empirical_factor=2.0):
    """Calculate the alpha param using the method Jasper, et al.

    :param n_eff: number of effective rotors
    :type n_eff: int
    :param zlj_dct: lennard-jones collision frequencies (cm-1?)
    :type zlj_dct: dict[float: float]
    :param empirical_factor: correction for using 1DME versus 2DM2
    :type empirical_factor: float
    """

    # Calculate the Lennard-Jones frequencies
    red_mass = (mass1 * mass2) / (mass1 + mass2)

    # Calculate Zalpha(Neff) at T = 300, 1000, 2000 K
    z_alphas_n_eff = _calculate_z_alpha_terms(n_eff, collider_set)

    # Calculate alpha = Zalpha(Neff) / Z(N) at T = 300, 1000, 2000 K
    # Empirical correction factor of (1/2) used for 1D Master Equations
    alpha_dct = {}
    for temp, z_alpha_n_eff in z_alphas_n_eff.items():
        zlj = troe_lj_collision_frequency(eps, sig, red_mass, temp)
        zlj *= 100**3  # conv. m^3/s to cm^3/ for below
        alpha_dct[temp] = (z_alpha_n_eff / zlj) / empirical_factor

    # Determine alpha and n for the e-down model
    edown_alpha, edown_n = _calculate_energy_down_exponent(alpha_dct)

    return edown_alpha, edown_n


def _calculate_z_alpha_terms(n_eff, collider_set):
    """Calculate the [Z*alpha](N_eff)"""

    def _z_alpha(n_eff, coeffs):
        """calculate an effective Z*alpha parameter"""
        return (
            coeffs[0] * n_eff ** (3)
            + coeffs[1] * n_eff ** (2)
            + coeffs[2] * n_eff ** (1)
            + coeffs[3]
        ) / 1.0e9

    # Read the proper coefficients from the moldriver dct
    coeff_dct = Z_ALPHA_EST_DCT.get(collider_set, None)
    if coeff_dct is not None:
        # Calculate the three alpha terms
        z_alpha_dct = {}
        for temp, coeffs in coeff_dct.items():
            z_alpha_dct[temp] = _z_alpha(n_eff, coeffs)

    return z_alpha_dct


def _calculate_energy_down_exponent(alpha_dct):
    """Calculate power n, for model: E_down = E_down_300 * (T/300)**n

    Does a least-squares for n to solve the linear equation
    ln(E_down/E_down_300) = [ln(T/300)] * n

    :param alpha_dct: temperature-dependent alpha parameters
    :type alpha_dct: dict[float: float]
    """

    assert 300 in alpha_dct, "Must have 300 K in alphas"

    # Set the edown alpha to the value at 300 K
    edown_alpha = alpha_dct[300]

    # Build vectors and matrices used for the fitting
    temps = numpy.array(list(alpha_dct.keys()), dtype=numpy.float64)
    alphas = numpy.array(list(alpha_dct.values()), dtype=numpy.float64)

    n_vec = numpy.log(temps / 300.0)
    coeff_mat = numpy.array([n_vec], dtype=numpy.float64)
    coeff_mat = coeff_mat.transpose()

    print("edown_vec test:", alphas, edown_alpha)

    edown_vec = numpy.log(alphas / edown_alpha)

    # Perform the least-squares fit
    theta = numpy.linalg.lstsq(coeff_mat, edown_vec, rcond=None)[0]

    # Set the the edown n value to the fitting parameter
    edown_n = theta[0]

    return edown_alpha, edown_n


# CALCULATE THE EFFECTIVE LENNARD-JONES SIGMA AND EPSILON
def lennard_jones_params(n_heavy, collider_set):
    """Returns in angstrom and cm-1.

    :param n_heavy: Number of heavy atoms for a species
    :type n_heavy: int
    """

    def _lj(n_heavy, coeff1, coeff2):
        """calculate Lennard-Jones parameter using estimation
        coefficents
        """
        return coeff1 * n_heavy ** (coeff2)

    # See if collider set is in dict with exact numbers
    # If not read the estimation coefficents from EST dct and calc. params
    params = LJ_DCT.get(collider_set)
    if params is not None:
        sig, eps = params
    else:
        coeffs = LJ_EST_DCT.get(collider_set, None)
        if coeffs is not None:
            # Calculate the effective sigma and epsilon values
            sig = _lj(n_heavy, coeffs[0], coeffs[1])
            eps = _lj(n_heavy, coeffs[2], coeffs[3])
        else:
            sig, eps = None, None

    # Convert the units to what they should be internally
    if sig is not None:
        sig *= phycon.ANG2BOHR
    if eps is not None:
        eps *= phycon.WAVEN2EH

    return sig, eps


# DETERMINE N_EFF USED FOR ALPHA AND LJ PARAM CALCULATIONS
def effective_rotor_count(geo):
    """Calculate an effective N parameter using the given parametrization.

    :param geo: geometry (Bohr)
    :type geo: automol geometry data structure
    :rtype: float
    """

    # Conver the geo to a graph
    gra = geom.graph(geo)
    symbs = geom.symbols(geo)

    # Count the rotors
    (
        n_pp,
        n_ps,
        n_pt,
        n_pq,
        n_ss,
        n_st,
        n_sq,
        n_tt,
        n_tq,
        n_qq,
        n_co,
        n_oo,
        n_ss_ring,
        n_rings,
    ) = _rotor_counts(gra, symbs)

    # print('    - Rotor Counts for N_eff:')
    # print('       N_pp:{}, N_ps:{}, N_pt:{}, N_pq:{}'.format(
    #     n_pp, n_ps, n_pt, n_pq))
    # print('       N_ss:{}, N_st:{}, N_sq:{}'.format(n_ss, n_st, n_sq))
    # print('       N_tt:{}, N_tq:{}'.format(n_tt, n_tq))
    # print('       N_qq:{}'.format(n_qq))
    # print('       N_co:{}, N_oo:{}'.format(n_co, n_oo))
    # print('       N_ss_ring:{}, N_rings:{}'.format(n_ss_ring, n_rings))

    # Use the rotor counts and the coefficients to calculate Neff
    c_pp_ps_ss, c_pt_st, c_pq_sq = 1.0, 2.0 / 3.0, 1.0 / 3.0
    c_tt_tq_qq, c_co_oo, c_ss_ring = 0.0, 1.0 / 3.0, 1.0 / 2.0
    n_eff = 1.0 + (
        c_pp_ps_ss * (n_pp + n_ps + n_ss)
        + c_pt_st * (n_pt + n_st)
        + c_pq_sq * (n_pq + n_sq)
        + c_tt_tq_qq * (n_tt + n_tq + n_qq)
        + c_co_oo * (n_co + n_oo)
        + c_ss_ring * n_ss_ring
        - n_rings
    )

    return n_eff


def _rotor_counts(gra, symbs):
    """Count up various types of bonds for a structure.

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
    rings = graph.rings(gra)
    ring_keys = set(itertools.chain(*graph.rings_atom_keys(gra)))
    n_rings = len(rings)

    # Loop over the bonds and count the number of atoms
    neighbors = graph.atoms_neighbor_atom_keys(gra)
    for bnd in graph.bond_keys(gra):
        key1, key2 = bnd
        spair = (symbs[key1], symbs[key2])
        if spair == ("C", "C"):
            # Figure out which neighbors are not hydrogen and count the number
            atom1_neighbors = neighbors[key1]
            numc1 = 0
            for neighbor1 in atom1_neighbors:
                if symbs[neighbor1] != "H":
                    numc1 += 1
            atom2_neighbors = neighbors[key2]
            numc2 = 0
            for neighbor2 in atom2_neighbors:
                if symbs[neighbor2] != "H":
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
        elif spair in (("C", "O"), ("O", "C")):
            n_co += 1
        elif spair == ("O", "O"):
            n_oo += 1

    # Compile counts into a tuple
    return (
        n_pp,
        n_ps,
        n_pt,
        n_pq,
        n_ss,
        n_st,
        n_sq,
        n_tt,
        n_tq,
        n_qq,
        n_co,
        n_oo,
        n_ss_ring,
        n_rings,
    )


# Rotational relaxation number
def rotational_relaxation_number(tgt_ich):
    """Get the rotational relaxation number at 298 K for a given species.
    Currently, this is read from internal library or set to 1.

    :param tgt_ich: InChI string of species
    :type tgt_ich: str
    :rtype: float
    """
    return ZROT_DCT.get(tgt_ich, 1.0)


# Determine which effective model series to use
def determine_collision_model_series(tgt_ich, bath_ich, collid_param):
    """For the collision between a given tgt and bath species, determine
    which effective series would be the most suitable model for
    determining the energy transfer parameters

    :param collid_param: select 'lj' or 'alpha' for parameter to assign
    """

    def _identify_effective_model(tgt_ich, bath_ich):
        """When values cannot be assigned for a collision, try and
        identify the best representative series to use for estimation
        """
        # Build the graph
        tgt_gra = geom.graph(chi.geometry(tgt_ich))

        # Identify the the target model
        if graph.is_radical_species(tgt_gra):
            # Determine if peroxy,hydroperoxy groups present to use RO2 series
            # otherwise just use alkyl radical series
            fgrp_cnt_dct = graph.functional_group_count_dct(tgt_gra)
            fgrps = set(fgrp for fgrp, count in fgrp_cnt_dct.items() if count > 0)
            _ro2_fgrps = {FunctionalGroup.PEROXY, FunctionalGroup.HYDROPEROXY}
            if _ro2_fgrps & fgrps:
                tgt_model = "peroxy"
            else:
                tgt_model = "1-alkyl"
        elif graph.is_hydrocarbon_species(tgt_gra):
            tgt_model = "n-alkane"
        else:
            # Set priority based on bond-dissociation energies
            # Loop through D0 dct (ordered by ene) and try to find func. grp
            tgt_model = None
            fgrp_cnt_dct = graph.functional_group_count_dct(tgt_gra)
            fgrps = set(fgrp for fgrp, count in fgrp_cnt_dct.items() if count > 0)
            for fgrp, model in D0_GRP_LST:
                if fgrp in fgrps:
                    tgt_model = model
                    break

            # Set target model to alkanes if nothing found
            if tgt_model is None:
                print("No info found for functional group {fgrps}")
                print("Using default model 'n-alkane'")
                tgt_model = "n-alkane"

        return frozenset({tgt_model, bath_ich})

    # First check if one should use standard numbers instead of estimating
    collider_set = frozenset({tgt_ich, bath_ich})

    # Go through a procedure to determine which model to use
    if collid_param == "lj":
        if collider_set not in LJ_DCT:
            # try to identify a model where possible
            collider_set = _identify_effective_model(tgt_ich, bath_ich)
    elif collid_param == "alpha":
        if {"InChI=1S/H2/h1H", "InChI=1S/H"} & collider_set:
            # Model cannot be set for these common situations
            # of H and H2 colliding with non Ar,N2 bath-gases
            collider_set = None
        else:
            # Last try to identify a model where possible
            collider_set = _identify_effective_model(tgt_ich, bath_ich)

    return collider_set
