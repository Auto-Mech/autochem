"""
  Use graph structures to identify chemical functional groups
"""

import itertools
from automol.graph._graph_base import atom_symbols
from automol.graph._graph_base import atom_symbol_idxs
from automol.graph._graph import rings
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import remove_atoms
from automol.graph._graph import full_isomorphism
from automol.graph._res import dominant_resonance
from automol.graph._res import resonance_dominant_radical_atom_keys
from automol.graph._res import bond_orders
from automol.graph._util import ring_idxs
from automol.graph._util import filter_idxs
from automol.graph._util import atom_idx_to_symb


class Fgroup():
    """ Functional groups
    """
    ALCOHOL = 'alcohol'
    PEROXY = 'peroxy'
    HYDROPEROXY = 'hydroperoxy'
    ETHER = 'ether'
    EPOXIDE = 'epoxide'
    ALDEHYDE = 'aldehyde'
    KETONE = 'ketone'
    ESTER = 'ester'
    CARBOX_ACID = 'carboxylic_acid'
    HALIDE = 'halide'
    THIOL = 'thiol'
    # TO ADD
    AMINE = 'amine'
    AMIDE = 'amide'
    NITRO = 'nitro'


def functional_group_dct(gra):
    """ Determine the functional groups for a given molecule.
    """

    # Build a dictionary by calling all the functional group functions
    # Certain smaller groups are removed when they are a part of larger groups
    peroxy_grps = peroxy_groups(gra)
    hydroperoxy_grps = hydroperoxy_groups(gra)
    ether_grps = ether_groups(gra)
    epoxide_grps = epoxide_groups(gra)
    carbox_acid_grps = carboxylic_acid_groups(gra)
    ester_grps = ester_groups(gra)
    ether_grps = ether_groups(gra, filterlst=ester_grps)
    alcohol_grps = alcohol_groups(gra, filterlst=carbox_acid_grps)
    aldehyde_grps = aldehyde_groups(gra, filterlst=carbox_acid_grps)
    ketone_grps = ketone_groups(gra, filterlst=carbox_acid_grps+ester_grps)
    # amine_grps = amine_groups(gra)
    amide_grps = amide_groups(gra)
    nitro_grps = nitro_groups(gra)
    halide_grps = halide_groups(gra)
    thiol_grps = thiol_groups(gra)

    # might have to filter it to remove ketone/oh if carbox acids are ther
    func_grp_dct = {
        Fgroup.PEROXY: peroxy_grps,
        Fgroup.HYDROPEROXY: hydroperoxy_grps,
        Fgroup.ETHER: ether_grps,
        Fgroup.EPOXIDE: epoxide_grps,
        Fgroup.CARBOX_ACID: carbox_acid_grps,
        Fgroup.ESTER: ester_grps,
        Fgroup.ALCOHOL: alcohol_grps,
        Fgroup.ALDEHYDE: aldehyde_grps,
        Fgroup.KETONE: ketone_grps,
        # Fgroup.AMINE: amine_grps,
        Fgroup.AMIDE: amide_grps,
        Fgroup.NITRO: nitro_grps,
        Fgroup.HALIDE: halide_grps,
        Fgroup.THIOL: thiol_grps
    }

    return func_grp_dct


# SEARCH FOR CERTAIN OVERARCHING MOLECULE TYPES
def hydrocarbon_species(gra):
    """ Determine if molecule is a hydrocarbon.

        :param uni_atoms: unique atomic symbols
        :type uni_atoms: tuple
        :rtype: boolean
    """
    return bool(set(_unique_atoms(gra)) <= {'C', 'H'})


def radical_species(gra):
    """ Determine if molecule is a radical species
    """
    return bool(resonance_dominant_radical_atom_keys(gra))


# SEARCH FOR REACTIVE SITES AND GROUPS
def alkene_sites(gra):
    """ Determine the location alkene groups
    """
    return bonds_of_type(gra, asymb1='C', asymb2='C', mbond=2)


def alkyne_sites(gra):
    """ Determine the location alkyne groups
    """
    return bonds_of_type(gra, asymb1='C', asymb2='C', mbond=3)


def alcohol_groups(gra, filterlst=()):
    """ Determine the location alcohol groups

        Returns a lsts of idxs of C-O-H groups
    """
    alc_grps = two_bond_idxs(gra, asymb1='C', cent='O', asymb2='H')
    alc_grps = filter_idxs(alc_grps, filterlst=filterlst)

    return alc_grps


def peroxy_groups(gra):
    """ Determine the location of hydroperoxy groups

        Returns a lsts of idxs of C-O-O groups
    """

    coo_r_grps = tuple()

    rad_idxs = resonance_dominant_radical_atom_keys(gra)

    coo_grps = two_bond_idxs(gra, asymb1='C', cent='O', asymb2='O')
    for coo_grp in coo_grps:
        c_idx, o1_idx, o2_idx = coo_grp
        if o2_idx in rad_idxs:
            coo_r_grps += ((c_idx, o1_idx, o2_idx),)

    return coo_r_grps


def hydroperoxy_groups(gra):
    """ Determine the location of hydroperoxy groups

        Returns a lsts of idxs of C-O-O-H groups
    """

    cooh_grps = tuple()

    coo_grps = two_bond_idxs(gra, asymb1='C', cent='O', asymb2='O')
    for coo_grp in coo_grps:
        c_idx, o1_idx, o2_idx = coo_grp
        o2_neighs = neighbors_of_type(gra, o2_idx, asymb='H')
        if o2_neighs:
            cooh_grps += ((c_idx, o1_idx, o2_idx, o2_neighs[0]),)

    return cooh_grps


def ether_groups(gra, filterlst=()):
    """ Determine the location of ether groups

        Returns a lsts of idxs of C-O-C groups
    """

    ether_grps = tuple()

    # Determing the indices of all rings in the molecule

    _ring_idxs = ring_idxs(rings(gra))

    coc_grps = two_bond_idxs(gra, asymb1='C', cent='O', asymb2='C')
    for coc_grp in coc_grps:
        c1_idx, o_idx, c2_idx = coc_grp
        if not _ring_idxs:
            ether_grps += ((c1_idx, o_idx, c2_idx),)
        else:
            for idxs in _ring_idxs:
                if not set(coc_grp) <= set(idxs):
                    ether_grps += ((c1_idx, o_idx, c2_idx),)

    ether_grps = filter_idxs(ether_grps, filterlst=filterlst)

    return ether_grps


def epoxide_groups(gra):
    """ Determine the location of epoxide groups

        Return C-O-C ring: only good for a 1,2-epoxide
    """

    epox_grps = tuple()

    # Determing the indices of all rings in the molecule
    _ring_idxs = ring_idxs(rings(gra))

    coc_grps = two_bond_idxs(gra, asymb1='C', cent='O', asymb2='C')
    for coc_grp in coc_grps:
        c1_idx, o_idx, c2_idx = coc_grp
        if _ring_idxs:
            for idxs in _ring_idxs:
                if set(coc_grp) <= set(idxs):
                    epox_grps += ((c1_idx, o_idx, c2_idx),)

    return epox_grps


def aldehyde_groups(gra, filterlst=()):
    """ Determine the location of aldehyde groups

        Returns C-O bond idxs
    """

    ald_grps = tuple()

    co_bonds = bonds_of_type(gra, asymb1='C', asymb2='O', mbond=2)
    for co_bond in co_bonds:
        c_idx, o_idx = co_bond
        c_neighs = neighbors_of_type(gra, c_idx, asymb='H')
        if c_neighs:
            ald_grps += ((c_idx, o_idx),)

    ald_grps = filter_idxs(ald_grps, filterlst=filterlst)

    return ald_grps


def ketone_groups(gra, filterlst=()):
    """ Determine the location of ketone groups

        Returns C-O bond idxs
    """

    ket_grps = tuple()

    co_bonds = bonds_of_type(gra, asymb1='C', asymb2='O', mbond=2)
    for co_bond in co_bonds:
        c_idx, o_idx = co_bond
        c_neighs = neighbors_of_type(gra, c_idx, asymb='H')
        if not c_neighs:
            ket_grps += ((c_idx, o_idx),)

    ket_grps = filter_idxs(ket_grps, filterlst=filterlst)

    return ket_grps


def ester_groups(gra):
    """ Determine the location of ester groups
        Likely identifies anhydrides as an ester
    """

    ester_grps = tuple()

    ether_grps = ether_groups(gra)
    ket_grps = ketone_groups(gra)

    for egrp, kgrp in itertools.product(ether_grps, ket_grps):
        c1_idx, o_idx, c2_idx = egrp

        ket_o_idx, ket_c_idx = None, None
        if c1_idx in kgrp:
            ket_o_idx = kgrp[1]
            ket_c_idx, term_c_idx = c1_idx, c2_idx
        elif c2_idx in kgrp:
            ket_o_idx = kgrp[1]
            ket_c_idx, term_c_idx = c2_idx, c1_idx

        if ket_o_idx is not None and ket_c_idx is not None:
            ester_grps += ((ket_o_idx, ket_c_idx, o_idx, term_c_idx),)

    return ester_grps


def carboxylic_acid_groups(gra):
    """ Determine the location of carboxylic acid groups
    """

    carbox_grps = tuple()

    alc_grps = alcohol_groups(gra)
    ket_grps = ketone_groups(gra)

    for agrp, kgrp in itertools.product(alc_grps, ket_grps):
        c_idx_set = list(set(agrp) & set(kgrp))
        if c_idx_set:
            ca_idx, oa_idx, ha_idx = agrp
            _, ok_idx = kgrp
            carbox_grps += ((ok_idx, ca_idx, oa_idx, ha_idx),)

    return carbox_grps


def amide_groups(gra):
    """ Determine the location of amide groups
    """

    amide_grps = tuple()

    ket_grps = ketone_groups(gra)
    noc_grps = two_bond_idxs(gra, asymb1='N', cent='O', asymb2='C')

    for noc_grp in noc_grps:
        n_idx, o_idx, c_idx = noc_grp
        if (c_idx, o_idx) in ket_grps:
            amide_grps += ((n_idx, o_idx, c_idx),)

    return amide_grps


def nitro_groups(gra):
    """ Determine the location alcohol groups

        Returns a lsts of idxs of NO2 groups
    """
    return two_bond_idxs(gra, asymb1='O', cent='N', asymb2='O')


def halide_groups(gra):
    """ Determine the location of halide groups
    """

    hal_grps = tuple()

    symb_idx_dct = atom_symbol_idxs(gra)

    for symb in ('F', 'Cl', 'Br', 'I'):
        hal_idxs = symb_idx_dct.get(symb, ())
        for hal_idx in hal_idxs:
            hal_neighs = neighbors_of_type(gra, hal_idx, asymb='C')
            hal_grps += ((hal_neighs[0], hal_idx),)

    return hal_grps


def thiol_groups(gra):
    """ Determine the location of thiol groups

        Returns a lsts of idxs of C-S-H groups
    """
    return two_bond_idxs(gra, asymb1='C', cent='S', asymb2='H')


# FIND GENERIC ATOM AND BOND GROUPS
def _unique_atoms(gra):
    """ Determine the symbols of unique atom types
    """
    symb_idx_dct = atom_symbol_idxs(gra)

    return tuple(symb_idx_dct.keys())


def chem_unique_atoms_of_type(gra, asymb):
    """ For the given atom type, determine the idxs of all the
         chemically unique atoms.
    """

    # Get the indices for the atom type
    symb_idx_dct = atom_symbol_idxs(gra)
    atom_idxs = symb_idx_dct[asymb]

    # Loop over each idx
    uni_idxs = tuple()
    uni_del_gras = []
    for idx in atom_idxs:

        # Remove the atom from the graph
        del_gra = remove_atoms(gra, [idx])

        # Test if the del_gra is isomorphic to any of the uni_del_gras
        new_uni = True
        for uni_del_gra in uni_del_gras:
            iso_dct = full_isomorphism(del_gra, uni_del_gra)
            if iso_dct:
                new_uni = False
                break

        # Add graph and idx to lst if del gra is unique
        if new_uni:
            uni_del_gras.append(del_gra)
            uni_idxs += (idx,)

    return uni_idxs


def bonds_of_type(gra, asymb1='C', asymb2='C', mbond=1):
    """ Determine the indices of all a specific bond
        specified by atom type and bond order
    """

    # Get the dict that relates atom indices to symbols
    idx_symb_dct = atom_symbols(gra)

    # Loop over all the bonds and build a list of ones that match
    _bonds_of_type = tuple()

    _bonds = bonds_of_order(gra, mbond=mbond)
    for bond in _bonds:
        idx1, idx2 = bond
        symb1, symb2 = idx_symb_dct[idx1], idx_symb_dct[idx2]
        if (symb1, symb2) in ((asymb1, asymb2), (asymb2, asymb1)):
            _bonds_of_type += ((idx1, idx2),)

    return _bonds_of_type


def bonds_of_order(gra, mbond=1):
    """ Determine the indices of a certain
    """

    # Build resonance graph to get the bond orders
    gra = dominant_resonance(gra)
    bond_order_dct = bond_orders(gra)

    mbond_idxs = tuple()
    for bond, order in bond_order_dct.items():
        if order == mbond:
            bnd1, bnd2 = bond
            mbond_idxs += ((bnd1, bnd2),)

    return mbond_idxs


def two_bond_idxs(gra, asymb1='H', cent='C', asymb2='H'):
    """ Determine triplet of idxs for describing bond

        idxs = (asymb1_idx, cent_idx, asymb2_idx)
    """

    grps = tuple()

    neigh_dct = atom_neighbor_keys(gra)
    idx_symb_dct = atom_symbols(gra)
    symb_idx_dct = atom_symbol_idxs(gra)

    cent_idxs = symb_idx_dct.get(cent, tuple())
    for cent_idx in cent_idxs:
        neighs = tuple(neigh_dct[cent_idx])
        neigh_symbs = atom_idx_to_symb(neighs, idx_symb_dct)
        if neigh_symbs == (asymb1, asymb2):
            grp_idxs = (neighs[0], cent_idx, neighs[1])
        elif neigh_symbs == (asymb2, asymb1):
            grp_idxs = (neighs[1], cent_idx, neighs[0])
        else:
            grp_idxs = ()

        if grp_idxs:
            grps += ((grp_idxs),)

    return grps


def neighbors_of_type(gra, aidx, asymb='H'):
    """ Get the neighbor indices for a certain type
    """

    idx_symb_dct = atom_symbols(gra)
    neighs = atom_neighbor_keys(gra)[aidx]
    neigh_symbs = atom_idx_to_symb(neighs, idx_symb_dct)

    idxs_of_type = tuple()
    for nidx, nsymb in zip(neighs, neigh_symbs):
        if nsymb == asymb:
            idxs_of_type += (nidx,)

    return idxs_of_type
