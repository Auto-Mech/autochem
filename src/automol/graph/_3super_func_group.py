"""
Assign species class based on CRECK classification
Implemented for aromatic species
Features: automatically recognize functional groups from species InChI
Identify them in graphs, and quantify how many there are
"""

import itertools

from .base import (
    # base functions
    implicit,
    # classification functions
    bonds_of_type,
    functional_group_dct,
)


def gra_has_grp(gra, grp):
    """filter species based on the presence of a functional group"""
    fc_grps_count = classify_species(gra)
    return bool(grp in fc_grps_count.keys())

def classify_species(gra):
    """uses the SuperFunctionalGroup 
    to classify a species according to species classes"""
    # call SuperFunctionalGroup
    fc_grps = SuperFunctionalGroup()
    fc_grps.assign_grps(gra)
    fc_grps.count_grps()
    return fc_grps.dct_count_grps
    # returns dct: number of groups for each type

RM_IF_OTHER_BASES = ["C3.DD-M", "C3.ST-M"]
#perhaps unnecessary since we filter other groups?
# TO ADD
# C3.D-M : C4H6, C3H6
# C3.DD-RSR : C3H3
# C4.DD-R : C4H5
# C4.DD-M : C4H6
# C4.DT-R : C4H3
# C3.D-RSR ; IC4H7, C4H71-3
# C5H4CCH2-RSR : FULVENALLENYL
# C5OH-RSR : C5H6OH
# A1,O2-M: BENZOQUINONE; FIND PROPER FUNCTIONAL GROUP
# A1CO-RSR: C6H5CO
# A1,H2-M : OC6H5CH3
# A1CH2O-R ? C6H5CH2O
# A1,C2H2-R; A1,C2H4-R => DIFFICILE PERCHé RAD NON è SU ANELLO
# A1,C3.DD-RSR: ELIMINA? reattività di A1CH2-RSR
BASE_GRP_DCT = {
    "C3.DD-M": "allene", # basis: should perhaps separate those without rings?
    "C3.ST-M": "propyne", 
    "C5-M": "cyclopentadiene",
    "C5O-M": "cyclopentadienone",
    "C5CH2-M": "fulvene",
    "FUR-M": "furan",
    "C5-RSR": "cyclopentadienyl",
    "C5H2-RSR": "cyclopentenyl",
    # SUBSTITUTED C5 RINGS
    "C5O-RSR": "cyclopentadienonyl",
    # AROMATICS
    "A1-M": "benzene",
    "A1-R": "phenyl",
    # SUBSTITUTED AROMATICS
    "A1CH2-RSR": "benzyl",
    # OXYGENATED AROMATICS
    "A1O-RSR": "phenoxy",
}
SUBSTITUENTS_GRP_DCT = {
    "OH": "alcohol",
    "CHO": "aldehyde",
    "CH3": "methyl",
    "C2H": "alkyne",
    "C2H3": "alkene",
    "C2H5": "alkane",
    "C3.DD": "allene",
    "C3.ST": "propyne",
    "OCH3": "alkoxy_oc",
    "C5": "cyclopentadiene",
    "A1": "aromatic", #as substituent: less restrictive
}

# POTENTIALLY, THE COMPOSITE GROUP LIST CAN BE MADE
# OF ALL THE STRUCTURES FROM THE BASE GROUP DICTIONARY
# COMBINED WITH ANY NUMBER AND TYPE OF SUBSTITUENTS.
# BUT THIS MAKES THE LIST SIMPLER AND MORE EFFECTIVE
# AND THE CODE FASTER
COMPOSITE_GRP_LIST = [
    # molecules - alkylated
    "C5,CH3-M",
    "A1,CH3-M",
    "A1,C2H-M",
    "A1,C2H3-M",
    "A1,C2H5-M",
    "A1,C3.DD-M",
    "A1,C3.ST-M",
    # to add: alkyl groups with radicals
    # molecules - oxygenated
    "C5,C5-M",
    "C5,OH-M",
    "A1,OH-M",
    "A1,OH,OH-M",
    "A1,OH,CH3-M",
    "A1,CH3,CH3-M",
    "A1,OH,CHO-M",
    "A1,OH,OCH3-M",
    "A1,CHO-M",
    "A1,OCH3-M",
    # radicals
    "C5,CH3-RSR",
    "C5,A1-RSR",
    "A1,CH3-R",
    "A1,OH-R",
    "A1O,OH-RSR",
    "A1O,CH3-RSR",
    "A1CH2,OH-RSR",
    "A1,C2H-R",
]


class SuperFunctionalGroup:
    """super functional groups composed of combinations of basic functional groups
    classification reflects that adopted in CRECK model for aromatic hydrocarbons
    """

    def __init__(
        self,
    ):
        self.sup_grps = {}
        self.grp_fct_dct = {}
        self.dct_count_grps = {}

    def assign_grps(self, gra):
        """assign sup groups to the graph provided
        """
        # call functional group dct
        self.grp_fct_dct = functional_group_dct(gra)

        # assign base groups
        for key, fct in BASE_GRP_DCT.items():
            self.sup_grps[key] = self.grp_fct_dct[fct]
        if len(self.sup_grps.keys()) > 1:
            for key in RM_IF_OTHER_BASES:
                self.sup_grps.pop(key, None)
        base_grps_0 = list(itertools.chain(*[
            grp for grp in self.sup_grps.values() if len(grp) > 0]))
        # if -RSR found but also -M and -R found: delete RSR
        # should have fixed this directly in functional group assignment
        # rtypes = [kk.split('-')[1] for kk in base_grps_0]
        # rsrtypes = [kk for kk in base_grps_0 if kk.split('-')[1] == 'RSR']
        # if len(rtypes) != len(rsrtypes): # 'RSR' but also other
        # e.g., 'R' and 'M' present: remove RSR types
        #     [self.sup_grps.pop(key, None) for key in rsrtypes]
        #     base_grps_0 = [base for base in base_grps_0 if '-RSR' not in base]
        # assign substituents
        subs_fct_dct = {}
        for key, fct in SUBSTITUENTS_GRP_DCT.items():
            subs_fct_dct[key] = self.grp_fct_dct[fct]
        # CHECK C6H5C2H2, C6H5C2H4!!
        # assign composite
        heavy_atms = list(implicit(gra)[0].keys())
        for comp_grp in COMPOSITE_GRP_LIST:
            base_and_subs, base_type = comp_grp.split("-")
            base, subs = (
                base_and_subs.split(",")[0] + "-" + base_type,
                base_and_subs.split(",")[1:],
            )

            base_grps = self.sup_grps[base]  # base groups to search substituents in
            for sub in subs:
                sub_grps = subs_fct_dct[sub]
                sub_grps_eff = ()

                # if the atoms of the substituent are part of (any) base group: skip
                for grp in sub_grps:
                    if not any(all(
                        atm in basei for atm in grp if atm in heavy_atms)
                        for basei in base_grps_0):
                        sub_grps_eff += (grp,)
                # intersection base+sub becomes the new base_grps;
                # filter by bond type, e.g., C-C, C-O..
                # with bonded_grps only: fails for OCH3
                # (CH2-O bonded to an aromatic would work too)
                base_grps = bonded_grps_checksymb(
                    gra, base_grps, sub_grps_eff, "C", sub[0])
            # add to dct
            self.sup_grps[comp_grp] = base_grps

    def count_grps(self):
        """count the assigned sup groups
        """
        self.dct_count_grps = {
            fgrp: len(grp_idx_lst)
            for fgrp, grp_idx_lst in self.sup_grps.items()
            if grp_idx_lst
        }


def bonded_grps(gra, grps1, grps2):
    """check if there is a bond between group1 and group2 of atoms in a graph
    return tuple of bonded groups
    grps1, grps2: tuple(tuple), groups of bonded atoms
    """
    heavy_atms = list(implicit(gra)[0].keys())
    grps = ()
    if len(grps1) > 0 and len(grps2) > 0:
        for grp1 in grps1:
            # keep only heavy atoms
            grp1 = tuple(atm for atm in grp1 if atm in heavy_atms)
            for grp2 in grps2:
                grp2 = tuple(
                    atm for atm in grp2 if atm in heavy_atms and atm not in grp1
                )
                possible_bonds = list(itertools.product(grp1, grp2))
                if any(frozenset(bond) in gra[1].keys() for bond in possible_bonds):
                    grp = grp1 + grp2
                    if sorted(grp) not in [sorted(grpi) for grpi in grps]:
                        grps += (grp,)

    return grps


def bonded_grps_checksymb(gra, grps1, grps2, symb1, symb2):
    """check if there is a bond between group1 and group2 of atoms in a graph
    return tuple of bonded groups
    grps1, grps2: tuple(tuple), groups of bonded atoms
    symb1, symb2: atom symbols of the bonded group sym1-sym2
    """
    heavy_atms = list(implicit(gra)[0].keys())
    correct_bonds = bonds_of_type(gra, symb1, symb2)
    grps = ()
    assigned_grps2 = () # make sure that each substituent is assigned to one group only
    if len(grps1) > 0 and len(grps2) > 0 and len(correct_bonds) > 0:
        for grp1 in grps1:
            # keep only heavy atoms
            grp1 = tuple(atm for atm in grp1 if atm in heavy_atms)
            for grp2 in grps2:
                grp2 = tuple(
                    atm for atm in grp2 if atm in heavy_atms and atm not in grp1
                )
                possible_bonds = list(itertools.product(grp1, grp2))
                possible_bonds += list(itertools.product(grp2, grp1))
                effective_bonds = list(
                    bond for bond in possible_bonds if frozenset(bond) in gra[1].keys()
                )
                if len(tuple(set(effective_bonds).intersection(correct_bonds))) > 0:
                    grp = grp1 + grp2
                    if (sorted(grp) not in [sorted(grpi) for grpi in grps]
                        and grp2 not in assigned_grps2):
                        grps += (grp,)
                        assigned_grps2 += (grp2,)
    return grps
