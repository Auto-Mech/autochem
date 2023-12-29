"""
  Use graph structures to identify chemical functional groups

  Note: Code requires explicit kekule graphs to work
"""

import itertools
from automol.graph.base._03kekule import kekule
from automol.graph.base._03kekule import radical_atom_keys
from automol.graph.base._02algo import branches
from automol.graph.base._02algo import rings_atom_keys
from automol.graph.base._02algo import isomorphism
from automol.graph.base._00core import atom_keys
from automol.graph.base._00core import bond_keys
from automol.graph.base._00core import atom_symbols
from automol.graph.base._00core import atom_symbol_keys
from automol.graph.base._00core import atoms_neighbor_atom_keys
from automol.graph.base._00core import bond_orders
from automol.graph.base._00core import remove_atoms
from automol.graph.base._00core import remove_bonds
from automol.graph.base._00core import ts_reactants_graph_without_stereo
from automol.graph.base._00core import subgraph
from automol.graph.base._00core import explicit
from automol.graph.base._11stereo import to_local_stereo
from automol.graph.base._11stereo import from_local_stereo


# # core functions
class FunctionalGroup():
    """ Functional groups
    """
    ALKENE = 'alkene'
    ALKOXY = 'alkoxy'
    ALCOHOL = 'alcohol'
    PEROXY = 'peroxy'
    HYDROPEROXY = 'hydroperoxy'
    ETHER = 'ether'
    CYCLIC_ETHER = 'cyclic_ether'
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
    METHYL = 'methyl'
    PHENYL = 'phenyl'


def functional_group_count_dct(gra):
    """ Return a dictionary that contains a count of the number
        of each of the functional groups in a species.

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: dict[str: int]
    """
    fgrp_dct = functional_group_dct(gra)
    return {fgrp: len(grp_idx_lst) for fgrp, grp_idx_lst in fgrp_dct.items()
            if grp_idx_lst}


def functional_group_dct(gra):
    """ Determine the functional groups for a given molecule.

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: dict[str: tuple(int)]
    """

    # Convert to explicit kekule graph for the functions to work
    gra = kekule(explicit(gra))

    # Build a dictionary by calling all the functional group functions
    # Certain smaller groups are removed when they are a part of larger groups
    alkene_grps = alkene_sites(gra)
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
    ketone_grps = ketone_groups(gra, filterlst=carbox_acid_grps+ester_grps)
    # amine_grps = amine_groups(gra)
    amide_grps = amide_groups(gra)
    nitro_grps = nitro_groups(gra)
    halide_grps = halide_groups(gra)
    thiol_grps = thiol_groups(gra)
    methyl_grps = methyl_groups(gra)
    phenyl_grps = phenyl_groups(gra)
    amine_grps = amine_groups(gra)
    # might have to filter it to remove ketone/oh if carbox acids are ther
    return {
        FunctionalGroup.ALKENE: alkene_grps,
        FunctionalGroup.ALKOXY: alkoxy_grps,
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
        FunctionalGroup.AMINE: amine_grps
    }


# # finders for overaching types
def is_hydrocarbon_species(gra):
    """ Determine if molecule is a hydrocarbon.

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: bool
    """
    return bool(set(_unique_atoms(gra)) <= {'C', 'H'})


def is_radical_species(gra):
    """ Determine if molecule is a radical species.

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: bool
    """
    return bool(radical_atom_keys(gra))


# # finders for reactive sites and groups
def alkene_sites(gra):
    """ Determine the location alkene groups

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: bool
    """
    alk_groups = ()
    cc2_bnds = bonds_of_type(gra, symb1='C', symb2='C', mbond=2)
    phenyl_grps = phenyl_groups(gra)
    for bnd in cc2_bnds:
        if not any(
                bnd in itertools.permutations(phen, r=2)
                for phen in phenyl_grps):
            alk_groups += (bnd,)
    return alk_groups


def alkyne_sites(gra):
    """ Determine the location alkyne groups

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: bool
    """
    return bonds_of_type(gra, symb1='C', symb2='C', mbond=3)


def alcohol_groups(gra, filterlst=()):
    """ Determine the location of alcohol groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-O-H atoms
        of the group: (C-idx, O-idx, H-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    alc_grps = two_bond_idxs(gra, symb1='C', cent='O', symb2='H')
    alc_grps = _filter_idxs(alc_grps, filterlst=filterlst)

    return alc_grps


def alkoxy_groups(gra):
    """ Determine the location of alkoxy groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-O atoms
        of the group: (C-idx, O-idx).

        Here the O-idx corresponds to a radical site.

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    alkox_grps = tuple()

    rad_keys = radical_atom_keys(gra)

    co_bonds = bonds_of_type(gra, symb1='C', symb2='O', mbond=1)
    for co_bond in co_bonds:
        c_idx, o_idx = co_bond
        if o_idx in rad_keys:
            alkox_grps += ((c_idx, o_idx),)

    return alkox_grps


def peroxy_groups(gra):
    """ Determine the location of peroxy groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-O-O atoms
        of the group: (C-idx, O-idx, O-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    coo_r_grps = tuple()

    rad_idxs = radical_atom_keys(gra)

    coo_grps = two_bond_idxs(gra, symb1='C', cent='O', symb2='O')
    for coo_grp in coo_grps:
        c_idx, o1_idx, o2_idx = coo_grp
        if o2_idx in rad_idxs:
            coo_r_grps += ((c_idx, o1_idx, o2_idx),)

    return coo_r_grps


def hydroperoxy_groups(gra):
    """ Determine the location of hydroperoxy groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-O-O-H atoms
        of the group: (C-idx, O-idx, O-idx, H-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    cooh_grps = tuple()

    coo_grps = two_bond_idxs(gra, symb1='C', cent='O', symb2='O')
    for coo_grp in coo_grps:
        c_idx, o1_idx, o2_idx, = coo_grp
        o2_neighs = neighbors_of_type(gra, o2_idx, symb='H')
        if o2_neighs:
            cooh_grps += ((c_idx, o1_idx, o2_idx, o2_neighs[0]),)

    return cooh_grps


def ether_groups(gra, filterlst=()):
    """ Determine the location of ether groups. The locations are
        specified as tuple of idxs indicating the C-O-C atoms
        of the group: (C-idx, O-idx, C-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    ether_grps = tuple()

    # Determing the indices of all rings in the molecule

    _ring_idxs = rings_atom_keys(gra)

    coc_grps = two_bond_idxs(gra, symb1='C', cent='O', symb2='C')
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
    """ Determine the location of cyclic ether groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-O-C atoms
        of the group: (C-idx, O-idx, C-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    cyc_ether_grps = tuple()

    # Determing the indices of all rings in the molecule
    _ring_idxs = rings_atom_keys(gra)

    coc_grps = two_bond_idxs(gra, symb1='C', cent='O', symb2='C')
    for coc_grp in coc_grps:
        c1_idx, o_idx, c2_idx = coc_grp
        if _ring_idxs:
            for idxs in _ring_idxs:
                if set(coc_grp) <= set(idxs):
                    cyc_ether_grps += ((c1_idx, o_idx, c2_idx),)

    return cyc_ether_grps


def aldehyde_groups(gra, filterlst=()):
    """ Determine the location of aldehyde groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-O atoms
        of the group: (C-idx, O-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    ald_grps = tuple()

    co_bonds = bonds_of_type(gra, symb1='C', symb2='O', mbond=2)
    for co_bond in co_bonds:
        c_idx, o_idx = co_bond
        c_neighs = neighbors_of_type(gra, c_idx, symb='H')
        if c_neighs:
            ald_grps += ((c_idx, o_idx),)

    ald_grps = _filter_idxs(ald_grps, filterlst=filterlst)

    return ald_grps


def ketone_groups(gra, filterlst=()):
    """ Determine the location of ketone groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-O atoms
        of the group: (C-idx, O-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    ket_grps = tuple()

    co_bonds = bonds_of_type(gra, symb1='C', symb2='O', mbond=2)
    for co_bond in co_bonds:
        c_idx, o_idx = co_bond
        c_neighs = neighbors_of_type(gra, c_idx, symb='H')
        if not c_neighs:
            ket_grps += ((c_idx, o_idx),)

    ket_grps = _filter_idxs(ket_grps, filterlst=filterlst)

    return ket_grps


def sulfanyl_groups(gra, filterlst=()):
    """ Determine the location of sulfanyl groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-O atoms
        of the group: (C-idx, O-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    sul_grps = tuple()

    so_bonds = bonds_of_type(gra, symb1='S', symb2='O', mbond=1)
    for so_bond in so_bonds:
        s_idx, o_idx = so_bond
        c_neighs = neighbors_of_type(gra, s_idx, symb='H')
        if not c_neighs:
            sul_grps += ((s_idx, o_idx),)

    sul_grps = _filter_idxs(sul_grps, filterlst=filterlst)

    return sul_grps


def ester_groups(gra):
    """ Determine the location of ester groups. The locations are
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
    """ Determine the location of ester groups. The locations are
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
    """ Determine the location of amide groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C(O)-O-H atoms
        of the group: (C-idx, O-idx, O-idx, H-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    amide_grps = tuple()

    ket_grps = ketone_groups(gra)
    noc_grps = two_bond_idxs(gra, symb1='N', cent='O', symb2='C')

    for noc_grp in noc_grps:
        n_idx, o_idx, c_idx = noc_grp
        if (c_idx, o_idx) in ket_grps:
            amide_grps += ((n_idx, o_idx, c_idx),)

    return amide_grps


def methyl_groups(gra):
    """ Determine the location of methyl groups. The locations are
        specified as tuple-of-tuple of idxs indicating the -CH3 atoms
        of the group: (C-idx, H-idx, H-idx, H-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """
    methyl_grps = tuple()
    ch_bonds = bonds_of_type(gra, symb1='C', symb2='H', mbond=1)
    for ch_x, ch_y, ch_z in itertools.combinations(ch_bonds, r=3):
        c_x, h_x = ch_x
        c_y, h_y = ch_y
        c_z, h_z = ch_z
        if c_x == c_y and c_x == c_z:
            methyl_grps += ((c_x, h_x, h_y, h_z),)
    return methyl_grps


def amine_groups(gra):
    """ Determine the location of amine groups. The locations are
        specified as tuple-of-tuple of idxs indicating the -CH3 atoms
        of the group: (N-idx, H-idx, H-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """
    amine_grps = tuple()
    nh_bonds = bonds_of_type(gra, symb1='N', symb2='H', mbond=1)
    for nh_x, nh_y in itertools.combinations(nh_bonds, r=2):
        n_x, h_x = nh_x
        n_y, h_y = nh_y
        if n_x == n_y:
            amine_grps += ((n_x, h_x, h_y),)
    return amine_grps


def nitro_groups(gra):
    """ Determine the location of nitro groups. The locations are
        specified as tuple-of-tuple of idxs indicating the O-N-O atoms
        of the group: (N-idx, O-idx, O-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    return two_bond_idxs(gra, symb1='O', cent='N', symb2='O',
                         allow_cent_conn=True)


def halide_groups(gra):
    """ Determine the location of halide groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-X atoms
        of the group: (C-idx, X-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """

    hal_grps = tuple()

    symb_idx_dct = atom_symbol_keys(gra)

    for symb in ('F', 'Cl', 'Br', 'I'):
        hal_idxs = symb_idx_dct.get(symb, ())
        for hal_idx in hal_idxs:
            hal_neighs = neighbors_of_type(gra, hal_idx, symb='C')
            if not hal_neighs:
                hal_neighs = neighbors_of_type(gra, hal_idx, symb='S')
                hal_neighs = neighbors_of_type(gra, hal_idx, symb='B')
                hal_neighs += neighbors_of_type(gra, hal_idx, symb='O')
                hal_neighs += neighbors_of_type(gra, hal_idx, symb='N')
            hal_grps += ((hal_neighs[0], hal_idx),)

    return hal_grps


def phenyl_groups(gra):
    """ Determine the location of phenyl groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-X atoms
        of the group: (C-idx, X-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """
    phenyl_grps = ()
    rngs_atm_keys = rings_atom_keys(gra)
    srt_rngs_atm_keys = [sorted(atm_keys) for atm_keys in rngs_atm_keys]
    cc2_grps = bonds_of_type(gra, symb1='C', symb2='C', mbond=2)
    for cc2_x, cc2_y, cc2_z in itertools.combinations(cc2_grps, r=3):
        if sorted(cc2_x + cc2_y + cc2_z) in srt_rngs_atm_keys:
            phenyl_grps += ((cc2_x + cc2_y + cc2_z),)
    return phenyl_grps


def thiol_groups(gra):
    """ Determine the location of thiol groups. The locations are
        specified as tuple-of-tuple of idxs indicating the C-S-H atoms
        of the group: (C-idx, S-idx, H-idx).

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(int)
    """
    return two_bond_idxs(gra, symb1='C', cent='S', symb2='H')


def ring_substituents(gra):
    """ Determine substituent groups on a ring
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
    func_grp_dct = functional_group_dct(gra)
    atm_symb_dct = atom_symbols(gra)
    ngb_atms_dct = atoms_neighbor_atom_keys(gra)
    for rng_keys in rngs_atm_keys:
        rngs_subst_gras[rng_keys] = {}
        for atm in rng_keys:
            groups = ()
            sub_atms = [
                natm for natm in ngb_atms_dct[atm] if natm not in rng_keys]
            for satm in sub_atms:
                if atm_symb_dct[satm] == 'H':
                    continue
                idented = False
                for grp, grp_idx_lst in func_grp_dct.items():
                    for idx_lst in grp_idx_lst:
                        if satm in idx_lst and atm:
                            if grp == FunctionalGroup.HALIDE:
                                groups += (atm_symb_dct[satm],)
                            else:
                                groups += (grp,)
                            idented = True
                            break
                    if idented:
                        break
                if not idented:
                    groups += (atm_symb_dct[satm] + '-chain',)
            rngs_subst_gras[rng_keys][atm] = groups
    return rngs_subst_gras


# # helper functions
def bonds_of_type(gra, symb1, symb2, mbond=1):
    """ Determine the indices of all a specific bond
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
    """ Determine the indices of all bonds in a molecule with
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
    """ Determine the triplet of indices of atoms of specified
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
    """ For a given atom, determine the indices of all the atoms
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
    """ Obtain the keys for atoms of the desired symbol that
        correspond to a radical site.

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :param aidx: index of atom for which to find neighbors
        :type aidx: int
        :param symb: symbols of desired atom types for neighbors
        :type symb: str
    """

    idx_symb_dct = atom_symbols(gra)

    rad_keys = ()
    for rad in radical_atom_keys(gra, sing_res=True):
        if idx_symb_dct[rad] == symb:
            rad_keys += (rad,)

    return rad_keys


def radical_dissociation_products(gra, pgra1):
    """ For a given species, determine the products of a dissociation
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


# helpers
def _unique_atoms(gra):
    """ Determine the symbols of unique atom types.

        :param gra: molecular graph
        :type gra: molecular graph data structure
        :rtype: tuple(str)
    """

    symb_idx_dct = atom_symbol_keys(gra)

    return tuple(symb_idx_dct.keys())


def _atom_idx_to_symb(idxs, idx_symb_dct):
    """ Convert a list of atom idxs (a1, a2, ..., an)
        to atom symbols
    """
    return tuple(idx_symb_dct[idx] for idx in idxs)


def _bond_idx_to_symb(idxs, idx_symb_dct):
    """ Convert a list of bond idxs ((a1, b1), (a2, b2), ..., (an, bn))
        to pairs of atom symbols
    """
    return tuple(
        (idx_symb_dct[idx1], idx_symb_dct[idx2]) for (idx1, idx2) in idxs
    )


def _filter_idxs(idxs_lst, filterlst=()):
    """ Filter out a tuple
    """

    filtered_lst = tuple()

    for idxs in idxs_lst:
        if not any(set(idxs) <= set(fidxs) for fidxs in filterlst):
            filtered_lst += (idxs,)

    return filtered_lst
