""" Handle symmetry factor stuff
"""
from . import geom, graph, reac, zmat
from .data import rotor
from .util import zmat_conv


# internal symmetry number
def symmetry_factors_from_sampling(geos, rotors, grxn=None):
    """Determines the internal symmetry number for a given conformer geometry
    by assessing a set of symmetrically similar structures that have been
    obtained by previous conformer sampling processes.

    (1) Explore saved conformers to find the list of similar conformers,
        i.e., those with a coulomb matrix and energy that are equivalent
        to those for the reference geometry.
    (2) Expand each of those similar conformers by applying
        rotational permutations to each of the terminal groups.
    (3) Count how many distinct distance matrices there are in
        the fully expanded conformer list.

    :param symm_geos: geometries symmetrically similar to one another
    """
    # First, get the end-group symmetry factor
    geo0 = geos[0]
    gra = geom.graph(geo0, stereo=False) if grxn is None else reac.ts_graph(grxn)
    end_sym_fac = end_group_symmetry_factor(geos[0], gra=gra)

    # Get indices for the torsions
    zma = rotor.rotors_zmatrix(rotors)
    zc_ = zmat.conversion_info(zma)
    tors_names = rotor.rotors_torsion_names(rotors, flat=True)
    tors_zidxs = [zmat.coordinate(zma, n) for n in tors_names]
    tors_gidxs = zmat_conv.relabel_zmatrix_key_sequence(zc_, tors_zidxs)

    # Ignore terminal hydrogens
    gra = graph.implicit(gra, atm_keys=graph.terminal_atom_keys(gra))
    idx_pool = graph.atom_keys(gra)
    tors_gidxs = [ixs for ixs in tors_gidxs if set(ixs) <= idx_pool]

    # For filtering out identical (not just symmetrically equivalent) structures
    def _are_equivalent(geo1, geo2):
        same_dist = geom.almost_equal_dist_matrix(
            geo1, geo2, thresh=3e-1, idxs=idx_pool
        )
        same_tors = geom.are_torsions_same(geo1, geo2, tors_gidxs)
        return same_dist and same_tors

    # Identify unique geometries
    uniq_geos = []
    for geo in geos:
        if not any(_are_equivalent(geo, g) for g in uniq_geos):
            uniq_geos.append(geo)

    int_sym_fac = len(uniq_geos) * end_sym_fac

    return int_sym_fac, end_sym_fac


def end_group_symmetry_factor(geo, gra=None) -> float:
    """Determine the symmetry factor for terminal groups in a geometry

    :param geo: A geometry
    :type geo: automol geom data structure
    :param gra: A graph describing connectivity, defaults to None
    :type gra: automol graph data structure, optional
    :returns: The symmetry factor for the terminal atoms
    :rtype: float
    """
    gra = geom.graph(geo, stereo=False) if gra is None else gra
    gra = graph.implicit(gra)

    term_nkey_dct = graph.terminal_atom_neighbors(gra, atom=False)

    bords_dct = graph.kekules_bond_orders_collated(gra)
    nhyd_dct = graph.atom_implicit_hydrogens(gra)

    factor = 1.0
    for term_key, term_nkey in term_nkey_dct.items():
        # If this can only be single-bonded, it is a rotational end-group For TS graphs,
        # this will exclude reacting bonds
        bkey = frozenset({term_key, term_nkey})
        is_single = set(bords_dct[bkey]) == {1}

        nhyd = nhyd_dct[term_key]
        if nhyd and is_single:
            factor *= nhyd

    return factor


def reduce_internal_symm(geo, int_symm, ext_symm, end_group_factor):
    """Reduce symm if external sym is 3??"""
    if ext_symm % 3 == 0 and end_group_factor > 1:
        if not graph.is_branched(geom.graph(geo)):
            int_symm = int_symm / 3

    return int_symm


def rotor_reduced_symm_factor(sym_factor, rotor_symms):
    """Decrease the overall molecular symmetry factor by the
    torsional mode symmetry numbers
    """
    for symm in rotor_symms:
        sym_factor /= symm

    return sym_factor


ICH_DCT = {"C": "InChI=1S/C", "O": "InChI=1S/O"}




def oxygenated_hydrocarbon_symm_num(
            geo, zrxn=None,
            account_for_enantiomer=True,
            radical_as_enantiomer=False):
    """determine the symmetry number of a CHNO molecule"
    :param geo: A geometry
    :type geo: automol geom data structure
    :param zrxn: reaction object
    :type zrxn: automol Reac object
    :param account_for_enantiomer: should symmetry number account for enantiomerism
    :type account_for_enantiomer: boolean
    :param radical_as_enantiomer: should radical site count as enantiomeric site?
    :type radical_as_enantiomer: boolean
    :returns: The internal and external symmetry factors
    :rtype: tuple(float, float)
    """
    def _atom_symmetry(atmi, atmj):
        """ Given axis atmi-atmj, returns the degree of symmetry
            on side atmi
        """
        atm_symm = 1
        # the if statement prevents double-doubling of symmetry
        # in cases of umbrella motions (see, symmetry_about_umbrella)
        if atmi not in set(rad_atms) ^ set(lone_pair_atms):
            ngb_pris = [
                priorities[atm] for atm in ngb_dct[atmi] if atm != atmj]
            if all(pri == ngb_pris[0] for pri in ngb_pris):
                atm_symm = len(ngb_pris)
        return atm_symm

    def _bond_symmetry(atmi, atmj):
        """ combined symmetry of ends of axis atmi-atmj
        """
        return max(
            _atom_symmetry(atmi, atmj),
            _atom_symmetry(atmj, atmi))

    def _is_enantiomeric():
        """ is enantiomeric given the given
            enantiomeric conditions
        """
        enantiomeric = False
        for (atm, ngbs) in ngb_dct.items():
            ngb_pris = [
                priorities[ngb] for ngb in ngbs]
            if not len(ngb_pris) == len(set(ngb_pris)):
                continue
            if len(ngb_pris) == 4:
                enantiomeric = True
            elif len(ngb_pris) == 3 and \
                 atm in rad_atms and \
                 radical_as_enantiomer:
                enantiomeric = True
            elif len(ngb_pris) == 3 and \
                 atm in lone_pair_atms:
                enantiomeric = True
        return enantiomeric

    def _symmetry_about_umbrella():
        """ the factor of symmetry about umbrella sites
        """
        umb_symm = 1.
        for umb_atm in set(rad_atms) ^ set(lone_pair_atms):
            if not any(umb_atm in rot_bnd for rot_bnd in rotat_bnds):
                continue
            ngb_pris = [
                priorities[atm] for atm in ngb_dct[umb_atm]]
            if len(ngb_pris) == 3 and \
               (len(ngb_pris) - len(set(ngb_pris))) == 1:
                umb_symm *= 2.
        return umb_symm

    # Get priorities and connectivity info from graph
    if zrxn is not None:
        gra = reac.ts_graph(zrxn)
    else:
        gra = geom.graph(geo)
    rotat_bnds = graph.rotational_bond_keys(gra)
    priorities = graph.canonical_priorities(gra)
    ngb_dct = graph.atoms_neighbor_atom_keys(gra)
    rad_atms = graph.radical_atom_keys(gra)
    lone_pair_atms = graph.lone_pair_atom_keys(gra)

    # Multiply the rotatable bond symmetries into the
    # internal symmetry number
    int_symm = 1.0
    for bnd in rotat_bnds:
        int_symm *= _bond_symmetry(*bnd)

    # Determine if the external symmetry number should
    # be half because of the enantiomer
    chiral_center = False
    if account_for_enantiomer:
        chiral_center = _is_enantiomeric()

    ext_symm = geom.external_symmetry_factor(
        geo, chiral_center=chiral_center)

    # Determine if the external symmetry number should
    # be doubled from symmetry around umbrella
    int_symm *= _symmetry_about_umbrella()

    return int_symm, ext_symm
