"""
Use graph structures to identify chemical functional groups

Note: Code requires explicit kekule graphs to work
"""

import itertools

from ._00core import (
    add_bonds,
    atom_keys,
    atom_neighbor_atom_keys,
    atom_symbol_keys,
    atom_symbols,
    atoms_bond_keys,
    atoms_neighbor_atom_keys,
    # atom_unpaired_electrons,
    bond_keys,
    bond_orders,
    explicit,
    # unsaturated_atom_keys
    from_data,
    implicit,
    remove_atoms,
    remove_bonds,
    subgraph,
    ts_reactants_graph_without_stereo,
    union,
)
from ._02algo import branches, isomorphism, rings_atom_keys
from ._03kekule import (
    atom_hybridizations_from_kekule,
    kekule,
    kekules,
    kekules_bond_orders_collated,
    # kekules_bond_orders,
    radical_atom_keys,
    radical_atom_keys_from_kekule,
)
from ._08canon import from_local_stereo, to_local_stereo


# # core functions
class FunctionalGroup:
    """Functional groups"""

    ALKANE = "alkane"
    ALKENE = "alkene"
    ALKYNE = "alkyne"
    ALLENE = "allene"
    PROPYNE = "propyne"
    ALLYL = "allyl"
    ALKOXY = "alkoxy"
    ALKOXY_OC = "alkoxy_oc"
    ALCOHOL = "alcohol"
    PEROXY = "peroxy"
    HYDROPEROXY = "hydroperoxy"
    ETHER = "ether"
    CYCLIC_ETHER = "cyclic_ether"
    ALDEHYDE = "aldehyde"
    KETONE = "ketone"
    ESTER = "ester"
    CARBOX_ACID = "carboxylic_acid"
    HALIDE = "halide"
    THIOL = "thiol"
    # TO ADD
    AMINE = "amine"
    AMIDE = "amide"
    NITRO = "nitro"
    METHYL = "methyl"
    PHENYL = "phenyl"
    # new
    BENZENE = "benzene"
    AROMATIC = "aromatic"
    BENZYL = "benzyl"
    FULV = "fulvene"
    CPDYL = "cyclopentadienyl"
    CPD = "cyclopentadiene"
    CPDONE = "cyclopentadienone"
    PHENOXY = "phenoxy"
    CPTYL = "cyclopentenyl"
    CPTOYL = "cyclopentadienonyl"
    FURAN = "furan"


def functional_group_count_dct(gra):
    """Return a dictionary that contains a count of the number
    of each of the functional groups in a species.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: dict[str: int]
    """
    fgrp_dct = functional_group_dct(gra)
    return {
        fgrp: len(grp_idx_lst) for fgrp, grp_idx_lst in fgrp_dct.items() if grp_idx_lst
    }


def functional_group_dct(gra):
    """Determine the functional groups for a given molecule.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: dict[str: tuple(int)]
    """

    # Convert to explicit kekule graph for the functions to work
    gra = kekule(explicit(gra))
    # rings to exclude
    c6_rings = ring_by_size_and_hyb(gra, hyb=2, ring_size=6, accept_notspX=2)
    # Build a dictionary by calling all the functional group functions
    # Certain smaller groups are removed when they are a part of larger groups
    aromatic_grps = aromatic_groups(gra)
    benzene_grps = benzene_groups(gra)
    phenyl_grps = phenyl_groups(gra)
    c5h5o_grps = c5h5o_groups(gra)
    cpd_grps = cyclopentadiene_groups(gra, filterlst=c5h5o_grps)
    fulv_grps = fulvene_groups(gra)
    cpdyl_grps = cyclopentadienyl_groups(gra)
    allyl_grps = allyl_groups_lowestspin(gra, filterlst=cpdyl_grps)
    cpdone_grps = cyclopentadienone_groups(gra)
    cptyl_grps = cyclopentenyl_groups(gra)
    furan_grps = furan_groups(gra)
    allene_grps = allene_sites(gra, filterlst=allyl_grps)
    propyne_grps = propyne_sites(gra)
    alkane_grps = alkane_sites(
        gra,
        filterlst=fulv_grps
        + cpdone_grps
        + cptyl_grps
        + c5h5o_grps
        + furan_grps
        + cpdyl_grps
        + c6_rings
        + allyl_grps
        + propyne_grps,
    )
    alkene_grps = alkene_sites(
        gra,
        filterlst=cpdyl_grps
        + cptyl_grps
        + c5h5o_grps
        + c6_rings
        + allyl_grps
        + allene_grps,
    )
    alkyne_grps = alkyne_sites(gra, filterlst=allyl_grps + allene_grps)
    alkoxy_grps = alkoxy_groups(gra)
    peroxy_grps = peroxy_groups(gra)
    hydroperoxy_grps = hydroperoxy_groups(gra)
    ether_grps = ether_groups(gra)
    cyc_ether_grps = cyclic_ether_groups(gra)
    carbox_acid_grps = carboxylic_acid_groups(gra)
    ester_grps = ester_groups(gra)
    ether_grps = ether_groups(gra, filterlst=ester_grps)
    alcohol_grps = alcohol_groups(gra, filterlst=carbox_acid_grps)
    aldehyde_grps = aldehyde_groups(gra, filterlst=carbox_acid_grps)
    ketone_grps = ketone_groups(gra, filterlst=carbox_acid_grps + ester_grps)
    alkoxy_oc_grps = alkoxy_OC_groups(gra, filterlst=aromatic_grps + alcohol_grps)
    # amine_grps = amine_groups(gra)
    amide_grps = amide_groups(gra)
    nitro_grps = nitro_groups(gra)
    halide_grps = halide_groups(gra)
    thiol_grps = thiol_groups(gra)
    methyl_grps = methyl_groups(gra, filterlst=propyne_grps)
    amine_grps = amine_groups(gra)
    benzyl_grps = benzyl_groups(gra, filterlst = benzene_grps + phenyl_grps)
    phenoxy_grps = phenoxy_groups(gra)

    # might have to filter it to remove ketone/oh if carbox acids are ther
    return {
        FunctionalGroup.ALKANE: alkane_grps,
        FunctionalGroup.ALKENE: alkene_grps,
        FunctionalGroup.ALKYNE: alkyne_grps,
        FunctionalGroup.ALLENE: allene_grps,
        FunctionalGroup.PROPYNE: propyne_grps,
        FunctionalGroup.ALLYL: allyl_grps,
        FunctionalGroup.ALKOXY: alkoxy_grps,
        FunctionalGroup.ALKOXY_OC: alkoxy_oc_grps,
        FunctionalGroup.PEROXY: peroxy_grps,
        FunctionalGroup.HYDROPEROXY: hydroperoxy_grps,
        FunctionalGroup.ETHER: ether_grps,
        FunctionalGroup.CYCLIC_ETHER: cyc_ether_grps,
        FunctionalGroup.CARBOX_ACID: carbox_acid_grps,
        FunctionalGroup.ESTER: ester_grps,
        FunctionalGroup.ALCOHOL: alcohol_grps,
        FunctionalGroup.ALDEHYDE: aldehyde_grps,
        FunctionalGroup.KETONE: ketone_grps,
        FunctionalGroup.AMIDE: amide_grps,
        FunctionalGroup.NITRO: nitro_grps,
        FunctionalGroup.HALIDE: halide_grps,
        FunctionalGroup.THIOL: thiol_grps,
        FunctionalGroup.METHYL: methyl_grps,
        FunctionalGroup.PHENYL: phenyl_grps,
        FunctionalGroup.AMINE: amine_grps,
        FunctionalGroup.BENZENE: benzene_grps,
        FunctionalGroup.AROMATIC: aromatic_grps,
        FunctionalGroup.BENZYL: benzyl_grps,
        FunctionalGroup.FULV: fulv_grps,
        FunctionalGroup.CPDYL: cpdyl_grps,
        FunctionalGroup.CPD: cpd_grps,
        FunctionalGroup.CPDONE: cpdone_grps,
        FunctionalGroup.PHENOXY: phenoxy_grps,
        FunctionalGroup.CPTYL: cptyl_grps,
        FunctionalGroup.CPTOYL: c5h5o_grps,
        FunctionalGroup.FURAN: furan_grps,
    }


# # finders for overaching types
def is_hydrocarbon_species(gra):
    """Determine if molecule is a hydrocarbon.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: bool
    """
    return bool(set(_unique_atoms(gra)) <= {"C", "H"})


def is_radical_species(gra):
    """Determine if molecule is a radical species.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: bool
    """
    return bool(radical_atom_keys(gra))


# # finders for reactive sites and groups
def alkane_sites(gra, filterlst=()):
    """Determine the location alkane groups
        as single C-C bonds (i.e., methane/methyl are excluded)
        exclude also: allyl grps, aromatic grps, any non-sp3 carbon
    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple
    """
    hyb_dct = atom_hybridizations_from_kekule(gra)
    non_sp3 = tuple(atm for atm in hyb_dct if hyb_dct[atm] != 3)
    cc1_bnds = bonds_of_type(gra, symb1="C", symb2="C", mbond=1)
    cc1_bnds = _filter_idxs(cc1_bnds, filterlst=filterlst)
    # delete bond if both are non-sp3 carbons
    # maybe classify according to : primary, secondary, tertiary?
    cc1_bnds = tuple(bnd for bnd in cc1_bnds if not all(atm in non_sp3 for atm in bnd))
    return cc1_bnds


def alkene_sites(gra, filterlst=()):
    """Determine the location alkene groups
        exclude also: allyl grps, aromatic grps
    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple
    """
    cc2_bnds = bonds_of_type(gra, symb1="C", symb2="C", mbond=2)
    cc2_bnds = _filter_idxs(cc2_bnds, filterlst=filterlst)
    return cc2_bnds


def alkyne_sites(gra, filterlst=()):
    """Determine the location alkyne groups
        exclude allyl resonances
    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple
    """
    cc3_bnds = bonds_of_type(gra, symb1="C", symb2="C", mbond=3)
    cc3_bnds = _filter_idxs(cc3_bnds, filterlst=filterlst)

    return cc3_bnds


def allene_sites(gra, filterlst=()):
    """determine the location of allene groups

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple
    """
    cc2_bnds = bonds_of_type(gra, symb1="C", symb2="C", mbond=2)
    double_bonds = _filter_idxs(cc2_bnds, filterlst=filterlst)

    return adjacent_grps(gra, double_bonds, double_bonds)


def propyne_sites(gra):
    """determine the location of propyne groups

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple
    """
    single_bonds = alkane_sites(gra)
    triple_bonds = alkyne_sites(gra)
    # check adjacent groups among those tested
    return adjacent_grps(gra, single_bonds, triple_bonds)


def alcohol_groups(gra, filterlst=()):
    """Determine the location of alcohol groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C-O-H atoms
    of the group: (C-idx, O-idx, H-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    alc_grps = two_bond_idxs(gra, symb1="C", cent="O", symb2="H")
    alc_grps = _filter_idxs(alc_grps, filterlst=filterlst)
    # return only the OH group
    alc_grps = tuple((atms[1], atms[2]) for atms in alc_grps)
    # tuple([(atms[1], atms[2]) for atms in alc_grps])
    return alc_grps


def alkoxy_groups(gra):
    """Determine the location of alkoxy groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C-O atoms
    of the group: (C-idx, O-idx).

    Here the O-idx corresponds to a radical site.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    alkox_grps = tuple()

    rad_keys = radical_atom_keys(gra)

    co_bonds = bonds_of_type(gra, symb1="C", symb2="O", mbond=1)
    for co_bond in co_bonds:
        c_idx, o_idx = co_bond
        if o_idx in rad_keys:
            alkox_grps += ((c_idx, o_idx),)

    return alkox_grps


def alkoxy_OC_groups(gra, filterlst=()):
    """Determine the location of alkoxy groups O-R. The locations are
    specified as tuple-of-tuple of idxs indicating the O-C atoms
    of the group: (O-idx, C-idx).
    the C in the bond must be sp3

    only molecular groups are accepted

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """
    alkox_grps = tuple()
    israd = is_radical_species(gra)
    hyb_dct = atom_hybridizations_from_kekule(gra)
    oc_bonds = bonds_of_type(gra, symb1="O", symb2="C", mbond=1)
    for oc_bond in oc_bonds:
        o_idx, c_idx = oc_bond
        if not israd and hyb_dct[c_idx] == 3:
            alkox_grps += ((o_idx, c_idx),)
    alkox_grps = _filter_idxs(alkox_grps, filterlst=filterlst)
    return alkox_grps


def peroxy_groups(gra):
    """Determine the location of peroxy groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C-O-O atoms
    of the group: (C-idx, O-idx, O-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    coo_r_grps = tuple()

    rad_idxs = radical_atom_keys(gra)

    coo_grps = two_bond_idxs(gra, symb1="C", cent="O", symb2="O")
    for coo_grp in coo_grps:
        c_idx, o1_idx, o2_idx = coo_grp
        if o2_idx in rad_idxs:
            coo_r_grps += ((c_idx, o1_idx, o2_idx),)

    return coo_r_grps


def hydroperoxy_groups(gra):
    """Determine the location of hydroperoxy groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C-O-O-H atoms
    of the group: (C-idx, O-idx, O-idx, H-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    cooh_grps = tuple()

    coo_grps = two_bond_idxs(gra, symb1="C", cent="O", symb2="O")
    for coo_grp in coo_grps:
        (
            c_idx,
            o1_idx,
            o2_idx,
        ) = coo_grp
        o2_neighs = neighbors_of_type(gra, o2_idx, symb="H")
        if o2_neighs:
            cooh_grps += ((c_idx, o1_idx, o2_idx, o2_neighs[0]),)

    return cooh_grps


def ether_groups(gra, filterlst=()):
    """Determine the location of ether groups. The locations are
    specified as tuple of idxs indicating the C-O-C atoms
    of the group: (C-idx, O-idx, C-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    ether_grps = tuple()

    # Determing the indices of all rings in the molecule

    _ring_idxs = rings_atom_keys(gra)

    coc_grps = two_bond_idxs(gra, symb1="C", cent="O", symb2="C")
    for coc_grp in coc_grps:
        c1_idx, o_idx, c2_idx = coc_grp
        if not _ring_idxs:
            ether_grps += ((c1_idx, o_idx, c2_idx),)
        else:
            for idxs in _ring_idxs:
                if not set(coc_grp) <= set(idxs):
                    ether_grps += ((c1_idx, o_idx, c2_idx),)

    ether_grps = _filter_idxs(ether_grps, filterlst=filterlst)

    return ether_grps


def cyclic_ether_groups(gra):
    """Determine the location of cyclic ether groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C-O-C atoms
    of the group: (C-idx, O-idx, C-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    cyc_ether_grps = tuple()

    # Determing the indices of all rings in the molecule
    _ring_idxs = rings_atom_keys(gra)

    coc_grps = two_bond_idxs(gra, symb1="C", cent="O", symb2="C")
    for coc_grp in coc_grps:
        c1_idx, o_idx, c2_idx = coc_grp
        if _ring_idxs:
            for idxs in _ring_idxs:
                if set(coc_grp) <= set(idxs):
                    cyc_ether_grps += ((c1_idx, o_idx, c2_idx),)

    return cyc_ether_grps


def aldehyde_groups(gra, filterlst=()):
    """Determine the location of aldehyde groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C-O atoms
    of the group: (C-idx, O-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    ald_grps = tuple()

    co_bonds = bonds_of_type(gra, symb1="C", symb2="O", mbond=2)
    for co_bond in co_bonds:
        c_idx, o_idx = co_bond
        c_neighs = neighbors_of_type(gra, c_idx, symb="H")
        if c_neighs:
            ald_grps += ((c_idx, o_idx),)

    ald_grps = _filter_idxs(ald_grps, filterlst=filterlst)

    return ald_grps


def ketone_groups(gra, filterlst=()):
    """Determine the location of ketone groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C-O atoms
    of the group: (C-idx, O-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    ket_grps = tuple()

    co_bonds = bonds_of_type(gra, symb1="C", symb2="O", mbond=2)
    for co_bond in co_bonds:
        c_idx, o_idx = co_bond
        c_neighs = neighbors_of_type(gra, c_idx, symb="H")
        if not c_neighs:
            ket_grps += ((c_idx, o_idx),)

    ket_grps = _filter_idxs(ket_grps, filterlst=filterlst)

    return ket_grps


def sulfanyl_groups(gra, filterlst=()):
    """Determine the location of sulfanyl groups. The locations are
    specified as tuple-of-tuple of idxs indicating the S-O atoms
    of the group: (S-idx, O-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    sul_grps = tuple()

    so_bonds = bonds_of_type(gra, symb1="S", symb2="O", mbond=1)
    for so_bond in so_bonds:
        s_idx, o_idx = so_bond
        c_neighs = neighbors_of_type(gra, s_idx, symb="H")
        if not c_neighs:
            sul_grps += ((s_idx, o_idx),)

    sul_grps = _filter_idxs(sul_grps, filterlst=filterlst)

    return sul_grps


def ester_groups(gra):
    """Determine the location of ester groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C(O)-O-C atoms
    of the group: (C-idx, O-idx, O-idx, C-idx).

    Likely identifies anhydrides as an ester.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    ester_grps = tuple()

    ether_grps = ether_groups(gra)
    ket_grps = ketone_groups(gra)

    for egrp, kgrp in itertools.product(ether_grps, ket_grps):
        c1_idx, o_idx, c2_idx = egrp

        ket_o_idx = ket_c_idx = term_c_idx = None
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
    """Determine the location of ester groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C(O)-O-H atoms
    of the group: (C-idx, O-idx, O-idx, H-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
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
    """Determine the location of amide groups. The locations are
    specified as tuple-of-tuple of idxs indicating the N-O-C atoms
    of the group: (N-idx, O-idx, C-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    amide_grps = tuple()

    ket_grps = ketone_groups(gra)
    noc_grps = two_bond_idxs(gra, symb1="N", cent="O", symb2="C")

    for noc_grp in noc_grps:
        n_idx, o_idx, c_idx = noc_grp
        if (c_idx, o_idx) in ket_grps:
            amide_grps += ((n_idx, o_idx, c_idx),)

    return amide_grps


def methyl_groups(gra, filterlst=()):
    """Determine the location of methyl groups. The locations are
    specified as tuple-of-tuple of idxs indicating the -CH3 atoms
    of the group: (C-idx, H-idx, H-idx, H-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """
    methyl_grps = tuple()
    ch_bonds = bonds_of_type(gra, symb1="C", symb2="H", mbond=1)
    for ch_x, ch_y, ch_z in itertools.combinations(ch_bonds, r=3):
        c_x, h_x = ch_x
        c_y, h_y = ch_y
        c_z, h_z = ch_z
        if c_x == c_y and c_x == c_z:
            methyl_grps += ((c_x, h_x, h_y, h_z),)

    methyl_grps = _filter_idxs(methyl_grps, filterlst=filterlst)

    return methyl_grps


def amine_groups(gra):
    """Determine the location of amine groups. The locations are
    specified as tuple-of-tuple of idxs indicating the -CH3 atoms
    of the group: (N-idx, H-idx, H-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """
    amine_grps = tuple()
    nh_bonds = bonds_of_type(gra, symb1="N", symb2="H", mbond=1)
    for nh_x, nh_y in itertools.combinations(nh_bonds, r=2):
        n_x, h_x = nh_x
        n_y, h_y = nh_y
        if n_x == n_y:
            amine_grps += ((n_x, h_x, h_y),)
    return amine_grps


def nitro_groups(gra):
    """Determine the location of nitro groups. The locations are
    specified as tuple-of-tuple of idxs indicating the O-N-O atoms
    of the group: (N-idx, O-idx, O-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    return two_bond_idxs(gra, symb1="O", cent="N", symb2="O", allow_cent_conn=True)


def halide_groups(gra):
    """Determine the location of halide groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C-X atoms
    of the group: (C-idx, X-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    hal_grps = tuple()

    symb_idx_dct = atom_symbol_keys(gra)

    for symb in ("F", "Cl", "Br", "I"):
        hal_idxs = symb_idx_dct.get(symb, ())
        for hal_idx in hal_idxs:
            hal_neighs = neighbors_of_type(gra, hal_idx, symb="C")
            if not hal_neighs:
                hal_neighs = neighbors_of_type(gra, hal_idx, symb="S")
                hal_neighs = neighbors_of_type(gra, hal_idx, symb="B")
                hal_neighs += neighbors_of_type(gra, hal_idx, symb="O")
                hal_neighs += neighbors_of_type(gra, hal_idx, symb="N")
            hal_grps += ((hal_neighs[0], hal_idx),)

    return hal_grps


def phenyl_groups(gra):
    """Determine the location of phenyl groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C-X atoms
    of the group: (C-idx, X-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """

    # A1-R: aromatic ring with sigma radical

    phenyl_grps = ()
    aro_grps = aromatic_groups(gra)
    if len(aro_grps) > 0:
        # check that radical is on atom of the group
        for aro_grp in aro_grps:
            # one of the carbons must have up to 2 neighbors / 2 bonds
            bd_keys = atoms_bond_keys(gra)
            if any(len(bd_keys[atm]) == 2 for atm in aro_grp):
                phenyl_grps += (aro_grp,)

    return phenyl_grps


def benzene_groups(gra):
    """determine location of benzene (molecular) groups as keus of heavy atoms

    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """
    benz_grps = ()
    aro_grps = aromatic_groups(gra)
    if len(aro_grps) > 0:
        # check if radical or molecule
        if not is_radical_species(gra):
            benz_grps = aro_grps

    return benz_grps


def benzyl_groups(gra, filterlst=()):
    """Determine location of benzyl groups as keys of heavy atoms
    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """
    bzyl_grps = tuple()
    if len(filterlst) > 0:
        return bzyl_grps
    # multiple kekules to check for resonance
    all_bd_ords = kekules_bond_orders_collated(gra)
    # get sp2 rings
    arom_grps = aromatic_groups(gra)
    cpdyl_grps = ring_by_size_and_hyb(gra, hyb=2, ring_size=5, accept_notspX=1)
    all_sp2_rings = tuple(itertools.chain(*arom_grps + cpdyl_grps))
    # check that the lateral group is a C, NOT part of an aromatic ring
    # AND that it can resonate on the bond with the ring
    # ngb_atms_dct = atoms_neighbor_atom_keys(gra)

    for aro_grp in arom_grps:
        for atm in aro_grp:
            sub_atm = [
                natm
                for natm in neighbors_of_type(gra, atm, "C")
                if natm not in all_sp2_rings
            ]
            if len(sub_atm) > 0:
                bd_ords = all_bd_ords[frozenset({atm, sub_atm[0]})]
                if 1 in bd_ords and 2 in bd_ords:
                    # can be of both order 1 and 2: resonant
                    bzyl_grps += (aro_grp + (sub_atm[0],),)
                    # add atom to the aromatic group
                    # # => this will define the benzyl group

    return bzyl_grps


def cyclopentenyl_groups(gra):
    """Determine location of cyclopentenyl-like groups as keys of heavy atoms
    e.g., CYC5H7
    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """
    cptyl_grps = tuple()
    # 5-membered rings with at least 1 double bond
    cpd_rings = ring_by_size_and_hyb(gra, hyb=2, ring_size=5, accept_notspX=3)
    # 5-membered rings with at least 3 sp3 carbons
    cpt_rings = ring_by_size_and_hyb(gra, hyb=3, ring_size=5, accept_notspX=2)
    # intersection
    cpd_rings_new = tuple(set(cpd_rings).intersection(cpt_rings))
    # allyl groups
    allyl_grps = allyl_groups(gra)
    # if ring has allyl-like resonance: it will be cptyl
    for ring in cpd_rings_new:
        # controlla che sia subset
        if any(set(allyl_grp).issubset(ring) for allyl_grp in allyl_grps):
            cptyl_grps += (ring,)

    return cptyl_grps


def c5h5o_groups(gra):
    """Determine location of c5h5o-like groups as keys of heavy atoms
    e.g., CYC5H5O InChI=1S/C5H5O/c6-5-3-1-2-4-5/h1-3H,4H2

    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """

    # 5-MEMBERED RINGS WITH AT LEAST ONE DOUBLE BOND
    # LATERAL GROUP IS OXYGEN AND IT NEEDS TO RESONATE WITH THE REST OF THE RING
    # OXYGEN WILL BE IDENTIFIED AS RADICAL IN THE REST OF THE STRUCTURE

    # multiple kekules to check for resonance
    all_bd_ords = kekules_bond_orders_collated(gra)
    # get sp2 rings
    cpd_rings = ring_by_size_and_hyb(gra, hyb=2, ring_size=5, accept_notspX=3)
    # check that the lateral group is a O
    # AND that it can resonate on the bond with the ring

    c5h5o_grps = resonant_oxygens(gra, cpd_rings, all_bd_ords)

    return c5h5o_grps


def cyclopentadienyl_groups(gra):
    """Determine location of cyclopentadienyl-like groups as keys of heavy atoms
    CHECK: C12H8, INDENYL, INDENE, C12H7
    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """
    cpdyl_grps = ()

    # get c5-memebered rings with at least two double bonds
    cpd_rings = ring_by_size_and_hyb(gra, hyb=2, ring_size=5, accept_notspX=1)

    # radicals- exclude those on atoms involved in double bonds
    # so you are automatically excluding also sigma radicals as in C12H7
    hyb_dct = atom_hybridizations_from_kekule(gra)
    all_rad_keys = radical_atom_keys(gra)
    pi_rad_keys = [rad for rad in all_rad_keys if hyb_dct[rad] != 2]
    # checks all implied in: c5 ring with at least 4 double bonds,
    # and check on radical that is not sigma
    # 1) resonance stabilization
    # 2) radical can be on the ring
    # 3) radical is not sigma

    for cpdyl_rng in cpd_rings:
        can_be_pi_rad = any(atm in pi_rad_keys for atm in cpdyl_rng)
        if can_be_pi_rad:
            cpdyl_grps += (cpdyl_rng,)

    return cpdyl_grps


def allyl_groups(gra):
    """Determine location of allyl groups
    if adjacent allyl-types, e.g., in benzyl., terbutyl
    check all resonant structures
    warning: in aromatic radicals, relatively unstable structures can also arise
    (e.g., in phenyl can also have an allylic-like form)
    if you want only the most stable resonance structure analyzed,
    use allyl_groups_lowestspin
    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys

    """

    def pairwise_bd_ords_check(bds_ords, ngb_rads):
        """check pairwise bonds of type single-double,
        and sigle has a res. radical on it"""
        flag = False
        for i, ord0 in enumerate(bds_ords[0]):
            bds_ords1 = bds_ords[1]
            ngb_rads0, ngb_rads1 = ngb_rads
            # one single and one double bond,
            # where the rad on the single bond is resonant
            if (
                ord0 == 1
                and ngb_rads0[i] == 1
                and bds_ords1[i] == 2
                and ngb_rads1[i] == 0
            ):
                flag = True
            elif (
                ord0 == 2
                and ngb_rads0[i] == 0
                and ngb_rads1[i] == 1
                and bds_ords1[i] == 1
            ):
                flag = True

        return flag

    all_rad_keys = radical_atom_keys(gra)
    # exclude radicals immediately
    if len(all_rad_keys) == 0:
        return ()

    allyl_grps = ()
    kgras = kekules(gra)
    kgra_rads = []
    for kgra in kgras:
        hyb_dct = atom_hybridizations_from_kekule(kgra)
        rad_dct = radical_atom_keys_from_kekule(kgra)
        # only sp3 carbons can have allyl-like resonance
        kgra_rads.append([rad for rad in rad_dct if hyb_dct[rad] == 3])
    all_bd_ords = kekules_bond_orders_collated(gra)

    carbons = atom_keys(gra, symb="C")
    # list of whether the carbon is a radical in each res. structure
    all_carb_israd = ()
    for carb in carbons:
        all_carb_israd += (
            tuple(1 if carb in kgra_rad else 0 for kgra_rad in kgra_rads),
        )

    for carb in carbons:
        # carbon neighbors
        ngbs = atom_neighbor_atom_keys(gra, carb, symb="C")
        if len(ngbs) < 2:
            continue

        res_allyl_atoms = (carb,)
        bds_ords = []
        ngb_rads = []
        for ngb in ngbs:
            bd_ords = all_bd_ords[frozenset({carb, ngb})]
            if 1 in bd_ords and 2 in bd_ords:
                res_allyl_atoms += (ngb,)
                bds_ords.append(bd_ords)
                ngb_rads.append(all_carb_israd[ngb])

        if len(res_allyl_atoms) == 3:
            flag = pairwise_bd_ords_check(bds_ords, ngb_rads)

        elif len(res_allyl_atoms) == 4:
            bds_ords_comb = list(itertools.combinations(bds_ords, 2))
            ngb_rads_comb = list(itertools.combinations(ngb_rads, 2))
            flags = []
            for bds_ords, ngb_rads in zip(bds_ords_comb, ngb_rads_comb):
                flags.append(pairwise_bd_ords_check(bds_ords, ngb_rads))
            flag = bool(True in flags)
        else:
            continue
        if flag is True:
            allyl_grps += (res_allyl_atoms,)

    return allyl_grps


def allyl_groups_lowestspin(gra, filterlst=()):
    """Determine location of allyl groups
    if adjacent allyl-types, e.g., in benzyl., terbutyl
    only for the most stable resonant structure
    warning: does not detect it for oxygenated structures such as phenoxy radical
    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys

    """
    allyl_grps = ()
    # hyb dct
    hyb_dct = atom_hybridizations_from_kekule(gra)
    # bond orders
    bd_ords = bond_orders(gra)
    # sp3 radical atom keys
    rad_dct = radical_atom_keys(gra)
    rad_dct = [rad for rad in rad_dct if hyb_dct[rad] == 3]

    # restrict to sp2 carbons in the middle of allyl group (sp2, not radical)
    sp2carbons = [
        carb
        for carb in atom_keys(gra, symb="C")
        if hyb_dct[carb] == 2 and carb not in rad_dct
    ]
    for carb in sp2carbons:
        # carbon neighbors
        ngbs = atom_neighbor_atom_keys(gra, carb, symb="C")
        if len(ngbs) < 2:
            continue

        res_allyl_atoms = (carb,)
        for ngb in ngbs:
            if (
                bd_ords[frozenset({carb, ngb})] == 2
                and ngb not in rad_dct
                or bd_ords[frozenset({carb, ngb})] == 1
                and ngb in rad_dct
            ):
                res_allyl_atoms += (ngb,)

        if len(res_allyl_atoms) >= 3:
            allyl_grps += (res_allyl_atoms,)

    allyl_grps = _filter_idxs(allyl_grps, filterlst=filterlst)
    return allyl_grps


def cyclopentadiene_groups(gra, filterlst=()):
    """Determine location of cyclopentadiene-like groups as keys of heavy atoms
    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """
    cpd_grps = tuple()

    if len(filterlst) > 0:
        return cpd_grps
    # get 5-memebered rings with at least two double bonds
    cpd_rings = ring_by_size_and_hyb(gra, hyb=2, ring_size=5, accept_notspX=1)
    # needs at least 1 sp3 atom
    cpt_rings = ring_by_size_and_hyb(gra, hyb=3, ring_size=5, accept_notspX=4)
    # intersection between these lists
    cpd_rings = tuple(set(cpd_rings).intersection(cpt_rings))
    # all atoms in ring have to be carbons
    atm_symb_dct = atom_symbols(gra)
    # exclude ring if it has radicals in it
    all_rad_keys = radical_atom_keys(gra)
    # exclude if lateral groups have only non-H atoms?
    for cpd_rng in cpd_rings:
        if not all(atm_symb_dct[atm] == "C" for atm in cpd_rng):
            continue
        if not any(atm in all_rad_keys for atm in cpd_rng):
            cpd_grps += (cpd_rng,)

    return cpd_grps


def furan_groups(gra):
    """Determine location of furan-like groups as keys of heavy atoms
    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """
    c5_grps = ()

    # get 5-memebered rings with at least two double bonds
    c5_rings = ring_by_size_and_hyb(gra, hyb=2, ring_size=5, accept_notspX=1)
    # all atoms in ring have to be carbons
    atm_symb_dct = atom_symbols(gra)
    # exclude ring if it has radicals in it
    all_rad_keys = radical_atom_keys(gra)
    # include if it has 4 carbons and one oxygen
    for c5_rng in c5_rings:
        symbs = [atm_symb_dct[atm] for atm in c5_rng]
        if symbs.count("C") != 4 or symbs.count("O") != 1:
            continue
        if not any(atm in all_rad_keys for atm in c5_rng):
            c5_grps += (c5_rng,)

    return c5_grps


def C5_dbl_bond(gra, side="CH2"):
    """searches for cyclopentadienone-like or fulvene-like structures
    side (str, optional): side group to search for. Defaults to 'CH2'.
    """
    c5_dbl_grps = ()

    # get 5-memebered rings with at least two double bonds
    # need all carbon atoms to be sp2
    cpd_rings = ring_by_size_and_hyb(gra, hyb=2, ring_size=5, accept_notspX=0)
    # all atoms in ring have to be carbons
    atm_symb_dct = atom_symbols(gra)
    # exclude ring if it has radicals in it
    all_rad_keys = radical_atom_keys(gra)
    for cpd_rng in cpd_rings:
        if not all(atm_symb_dct[atm] == "C" for atm in cpd_rng):
            continue
        if any(atm in all_rad_keys for atm in cpd_rng):
            continue
        for atm_rng in cpd_rng:
            sub_atm = []
            if side == "O":
                # check that first atom of lateral substituent is an oxygen
                sub_atm = neighbors_of_type(gra, atm_rng, "O")
            elif side == "CH2":
                sub_atm = [
                    natm
                    for natm in neighbors_of_type(gra, atm_rng, "C")
                    if len(neighbors_of_type(gra, natm, "H")) == 2
                ]
            if len(sub_atm) > 0:
                c5_dbl_grps += (cpd_rng + (sub_atm[0],),)
    return c5_dbl_grps


def cyclopentadienone_groups(gra):
    """Determine location of cyclopentadienone-like groups as keys of heavy atoms
       cpd rings with all sp2 carbons and first atom of one lateral group is an oxygen
    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """
    cpdo_grps = C5_dbl_bond(gra, side="O")

    return cpdo_grps


def fulvene_groups(gra):
    """Determine location of fulvalene-like groups as keys of heavy atoms
       cpd rings with all sp2 carbons and first atom of
       one lateral group is an sp2 carbon

    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """
    flvl_grps = C5_dbl_bond(gra, side="CH2")

    return flvl_grps


def phenoxy_groups(gra):
    """Determine location of phenoxy groups as keys of heavy atoms
    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """
    # multiple kekules to check for resonance
    all_bd_ords = kekules_bond_orders_collated(gra)
    # get sp2 rings
    arom_grps = aromatic_groups(gra)
    # check that the lateral group is a O
    # AND that it can resonate on the bond with the ring
    phoxy_grps = resonant_oxygens(gra, arom_grps, all_bd_ords)

    return phoxy_grps


def aromatic_groups(gra):
    """Determine location of aromatic groups as keys of heavy atoms
    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    rtype: tuple(tuple) with non-H keys
    """
    arom_grps = ring_by_size_and_hyb(gra, hyb=2, ring_size=6)

    return arom_grps


def thiol_groups(gra):
    """Determine the location of thiol groups. The locations are
    specified as tuple-of-tuple of idxs indicating the C-S-H atoms
    of the group: (C-idx, S-idx, H-idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """
    return two_bond_idxs(gra, symb1="C", cent="S", symb2="H")


def ring_by_size_and_hyb(gra, hyb=2, ring_size=6, accept_notspX=0):
    """Determine location of rings with atoms with h
        hybridization hyb, of size ring_size, with up to accept_notspX
        atoms with "wrong" hybridization

    Args:
    param gra: molecular graph (kekule)
    type gra: tuple(dct)
    param hyb: desired hybridization of the atoms in the ring
    type hyb: int
    param ring_size: N of atoms in the ring
    type ring_size: int
    param accept_notspX: how many atoms with undesired hybridization are accepted?
    type accept_notspX: int
    rtype: tuple(tuple) with non-H keys
    """
    rng_keys_lst = rings_atom_keys(gra)
    hyb_dct = atom_hybridizations_from_kekule(gra)
    ring_grps = ()
    for rng_keys in rng_keys_lst:
        if len(rng_keys) == ring_size and sum(hyb_dct[k] == hyb for k in rng_keys) >= (
            ring_size - accept_notspX
        ):
            ring_grps += (rng_keys,)

    return ring_grps


def ring_substituents(gra):
    """Determine substituent groups on a ring
    to produce a graph of graphs where the top level
    key of a ring_gra is the order of the atm keys
    that define the ring
    aka (a1, a2, a3, a4, a5, a6) a1 is the 0th position of the ring
    so a3-a5 have a 1-3 interaction.  The nested dictionary has
    atm key: tuple of groups

    (a1, a2, a3, a4, a5, a6): {
        a1: (group1, )
        a2: (group1, group2, )
        a3: ()
        a4: ()
        a5: ()
        a6: (group1,)
    }
    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(int)
    """
    rngs_subst_gras = {}
    rngs_atm_keys = rings_atom_keys(gra)
    func_grp_dct = functional_group_dct(gra)  # find all possible functional groups
    atm_symb_dct = atom_symbols(gra)
    ngb_atms_dct = atoms_neighbor_atom_keys(gra)
    for rng_keys in rngs_atm_keys:
        rngs_subst_gras[rng_keys] = {}
        for atm in rng_keys:
            groups = ()
            sub_atms = [natm for natm in ngb_atms_dct[atm] if natm not in rng_keys]
            for satm in sub_atms:
                if atm_symb_dct[satm] == "H":
                    continue
                idented = False
                for grp, grp_idx_lst in func_grp_dct.items():
                    for idx_lst in grp_idx_lst:
                        if (
                            satm in idx_lst and atm
                        ):  # if substituent is in the identified group: add
                            if grp == FunctionalGroup.HALIDE:
                                groups += (atm_symb_dct[satm],)
                            else:
                                groups += (grp,)
                            idented = True
                            break
                    if idented:
                        break
                if not idented:
                    groups += (atm_symb_dct[satm] + "-chain",)
            rngs_subst_gras[rng_keys][atm] = groups
    return rngs_subst_gras


# # helper functions
def adjacent_grps(gra, grps1, grps2):
    """check if grp1 and grp2 are adjacent
    by detecting an atom in common
    return grps of bonded atoms without redundancy
    """
    heavy_atms = list(implicit(gra)[0].keys())  # keep only heavy atoms
    grps = ()
    if len(grps1) > 0 and len(grps2) > 0:
        for grp1 in grps1:
            grp1 = tuple(atm for atm in grp1 if atm in heavy_atms)
            for grp2 in grps2:
                grp2 = tuple(atm for atm in grp2 if atm in heavy_atms)
                if len(tuple(set(grp1).intersection(grp2))) > 0 and grp1 != grp2:
                    grpadd = tuple(set(grp1 + grp2))
                    if grpadd not in grps:
                        grps += (grpadd,)

    return grps


def resonant_oxygens(gra, grps, all_bd_ords):
    """find lateral oxygens resonant with the grps
    gra: species graph
    grps: groups to be checked for oxygen resonance
    all_bd_ords: bond orders for all kekule structures
    """
    res_grps = ()
    for grp in grps:
        for atm in grp:
            sub_atm = list(neighbors_of_type(gra, atm, "O"))
            if len(sub_atm) > 0:
                bd_ords = all_bd_ords[frozenset({atm, sub_atm[0]})]
                if 1 in bd_ords and 2 in bd_ords:
                    # can be of both order 1 and 2: resonant
                    res_grps += (grp + (sub_atm[0],),)
                    # add atom to the aromatic group
                    # # => this will define the phenoxy group
    return res_grps


def bonds_of_type(gra, symb1, symb2, mbond=1):
    """Determine the indices of all a specific bond
    specified by atom type and bond order.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :param symb1: symbol of atom 1 in the bond
    :type symb1: str
    :param symb2: symbol of atom 2 in the bond
    :type symb2: str
    :param mbond: bond order of desired bond type
    :type mbond: int
    :rtype: tuple(int)
    """

    # Get the dict that relates atom indices to symbols
    idx_symb_dct = atom_symbols(gra)

    # Loop over all the bonds and build a list of ones that match
    _bonds_of_type = tuple()

    _bonds = bonds_of_order(gra, mbond=mbond)
    for bond in _bonds:
        idx1, idx2 = bond
        _symb1, _symb2 = idx_symb_dct[idx1], idx_symb_dct[idx2]
        if (_symb1, _symb2) == (symb1, symb2):
            _bonds_of_type += ((idx1, idx2),)
        elif (_symb1, _symb2) == (symb2, symb1):
            _bonds_of_type += ((idx2, idx1),)

    return _bonds_of_type


def bonds_of_order(gra, mbond=1):
    """Determine the indices of all bonds in a molecule with
    the specified bond order.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :param mbond: bond order of desired bond type
    :type mbond: int
    """

    bond_order_dct = bond_orders(gra)

    mbond_idxs = tuple()
    for bond, order in bond_order_dct.items():
        if order == mbond:
            bnd1, bnd2 = bond
            mbond_idxs += ((bnd1, bnd2),)

    return mbond_idxs


def two_bond_idxs(gra, symb1, cent, symb2, allow_cent_conn=False):
    """Determine the triplet of indices of atoms of specified
    types that are connected in a chain by two bonds:
    (symb1_idx, cent_idx, symb2_idx).

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :param symb1: symbol of atom at one end of chain
    :type symb1: str
    :param cent: symbol of atom in the middle of a chain
    :type cent: str
    :param symb2: symbol of atom at other end of chain
    :type symb2: str
    """

    grps = tuple()

    neigh_dct = atoms_neighbor_atom_keys(gra)
    idx_symb_dct = atom_symbols(gra)
    symb_idx_dct = atom_symbol_keys(gra)

    cent_idxs = symb_idx_dct.get(cent, tuple())
    for cent_idx in cent_idxs:
        grp_idxs = None
        neighs = tuple(neigh_dct[cent_idx])
        neigh_symbs = _atom_idx_to_symb(neighs, idx_symb_dct)
        if neigh_symbs == (symb1, symb2):
            grp_idxs = (neighs[0], cent_idx, neighs[1])
        elif neigh_symbs == (symb2, symb1):
            grp_idxs = (neighs[1], cent_idx, neighs[0])
        elif allow_cent_conn and symb1 in neigh_symbs and symb2 in neigh_symbs:
            idx_1 = None
            idx_2 = None
            for idx in neighs:
                if idx_symb_dct[idx] == symb1 and idx_1 is None:
                    idx_1 = idx
                elif idx_symb_dct[idx] == symb2:
                    idx_2 = idx
            if idx_1 is not None and idx_2 is not None:
                grp_idxs = (idx_1, cent_idx, idx_2)
        else:
            grp_idxs = None

        if grp_idxs is not None:
            grps += ((grp_idxs),)

    return grps


def neighbors_of_type(gra, aidx, symb):
    """For a given atom, determine the indices of all the atoms
    which neighbor it that are of the type specified.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :param aidx: index of atom for which to find neighbors
    :type aidx: int
    :param symb: symbols of desired atom types for neighbors
    :type symb: str
    """

    idx_symb_dct = atom_symbols(gra)
    neighs = atoms_neighbor_atom_keys(gra)[aidx]
    neigh_symbs = _atom_idx_to_symb(neighs, idx_symb_dct)

    idxs_of_type = tuple()
    for nidx, nsymb in zip(neighs, neigh_symbs):
        if nsymb == symb:
            idxs_of_type += (nidx,)

    return idxs_of_type


def radicals_of_type(gra, symb):
    """Obtain the keys for atoms of the desired symbol that
    correspond to a radical site.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :param symb: symbols of desired radical sites
    :type symb: str
    """

    idx_symb_dct = atom_symbols(gra)

    rad_keys = ()
    for rad in radical_atom_keys(gra, sing_res=True):
        if idx_symb_dct[rad] == symb:
            rad_keys += (rad,)

    return rad_keys


def radical_dissociation_products(gra, pgra1):
    """For a given species, determine the products of a dissociation
    occuring around a radical site. We assume one of the
    dissociation products is known, and we attempt to find the
    corresponding product.

    Currently, we assume that the input pgra1 is appropriately
    stereolabeled.

    :param gra: species undergoing dissociation
    :type gra: automol.graph object
    :param pgra1: one of the known products of dissociation
    :type pgra1: automol.graph object
    :rtype: tuple(automol.graph.object)
    """

    gra = ts_reactants_graph_without_stereo(gra)

    # Attempt to find a graph of product corresponding to pgra1
    pgra2 = None
    for rad in radical_atom_keys(gra, sing_res=True):
        for adj in atoms_neighbor_atom_keys(gra)[rad]:
            for group in branches(gra, adj, stereo=False):
                if isomorphism(group, pgra1, backbone_only=True):
                    pgra2 = remove_atoms(gra, atom_keys(group), stereo=False)
                    if bond_keys(group) in pgra2:
                        pgra2 = remove_bonds(pgra2, bond_keys(group))

    # If pgra2 is ID'd, rebuild the two product graphs with stereo labels
    if pgra2 is not None:
        keys2 = atom_keys(pgra2)
        idx_gra = to_local_stereo(gra)
        idx_pgra2 = subgraph(idx_gra, keys2, stereo=True)
        pgra2 = from_local_stereo(idx_pgra2)

    return pgra1, pgra2


def possible_radical_dissociation_sources(prd_gras):
    """Given product graphs, determine possible radical dissociation sources that could
    have produced them, if any.

    :param gras: Product graphs
    :return: Product graphs, with combined unstable product, if any
    """
    oh_gra0 = from_data({0: "O", 1: "H"}, [(0, 1)])

    prd_gras0 = []
    oh_gra = None
    oh_key = None
    for gra in prd_gras:
        map_dct = isomorphism(oh_gra0, gra, backbone_only=True)
        if map_dct and not oh_gra:
            oh_gra = gra
            oh_key = map_dct.get(0, None)
        else:
            prd_gras0.append(gra)

    prds_gras = []
    for idx, src_gra in enumerate(prd_gras0):
        prd_gras = prd_gras0.copy()
        co_bonds = bonds_of_type(kekule(src_gra), symb1="C", symb2="O", mbond=2)
        for _, co_key in co_bonds:
            oo_bkey = frozenset({oh_key, co_key})
            src_gra = union(src_gra, oh_gra)
            src_gra = add_bonds(src_gra, [oo_bkey], ord_dct={oo_bkey: 0.8})
            prd_gras[idx] = src_gra
            prds_gras.append(prd_gras)

    return prds_gras


# helpers
def _unique_atoms(gra):
    """Determine the symbols of unique atom types.

    :param gra: molecular graph
    :type gra: molecular graph data structure
    :rtype: tuple(str)
    """

    symb_idx_dct = atom_symbol_keys(gra)

    return tuple(symb_idx_dct.keys())


def _atom_idx_to_symb(idxs, idx_symb_dct):
    """Convert a list of atom idxs (a1, a2, ..., an)
    to atom symbols
    """
    return tuple(idx_symb_dct[idx] for idx in idxs)


def _bond_idx_to_symb(idxs, idx_symb_dct):
    """Convert a list of bond idxs ((a1, b1), (a2, b2), ..., (an, bn))
    to pairs of atom symbols
    """
    return tuple((idx_symb_dct[idx1], idx_symb_dct[idx2]) for (idx1, idx2) in idxs)


def _filter_idxs(idxs_lst, filterlst=()):
    """Filter out a tuple"""
    filterlst = tuple(itertools.chain(*filterlst))
    filtered_lst = tuple(
        bnd for bnd in idxs_lst if not any(atm in filterlst for atm in bnd)
    )
    return filtered_lst


def _filter_idxs_old(idxs_lst, filterlst=()):
    """Filter out a tuple"""

    filtered_lst = tuple()

    for idxs in idxs_lst:
        if not any(set(idxs) <= set(fidxs) for fidxs in filterlst):
            filtered_lst += (idxs,)

    return filtered_lst
