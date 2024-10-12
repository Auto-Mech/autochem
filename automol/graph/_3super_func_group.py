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


BASE_GRP_DCT = {
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
    "A1,C3.DD-M",
    "A1,C3.ST-M",
    # molecules - oxygenated
    "C5,OH-M",
    "A1,OH-M",
    "A1,OH,OH-M",
    "A1,OH,CHO-M",
    "A1,OH,OCH3-M",
    "A1,CHO-M",
    "A1,OCH3-M",
    # radicals
    "C5,CH3-RSR",
    "A1,CH3-R",
    "A1,OH-R",
    "A1O,OH-RSR",
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

        # assign substituents
        subs_fct_dct = {}
        for key, fct in SUBSTITUENTS_GRP_DCT.items():
            subs_fct_dct[key] = self.grp_fct_dct[fct]

        # CH3CK C6H5C2H2, C6H5C2H4!!
        # assign composite
        for comp_grp in COMPOSITE_GRP_LIST:
            base_and_subs, base_type = comp_grp.split("-")
            base, subs = (
                base_and_subs.split(",")[0] + "-" + base_type,
                base_and_subs.split(",")[1:],
            )
            base_grps = self.sup_grps[base]  # base groups to search substituents in
            for sub in subs:
                sub_grps = subs_fct_dct[sub]
                # intersection becomes the new base_grps;
                # filter by bond type, e.g., C-C, C-O..
                # with bonded_grps only: fails for OCH3
                # (CH2-O bonded to an aromatic would work too)
                base_grps = bonded_grps_checksymb(gra, base_grps, sub_grps, "C", sub[0])
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
    if len(grps1) > 0 and len(grps2) > 0 and len(correct_bonds) > 0:
        for grp1 in grps1:
            # keep only heavy atoms
            grp1 = tuple(atm for atm in grp1 if atm in heavy_atms)
            for grp2 in grps2:
                grp2 = tuple(
                    atm for atm in grp2 if atm in heavy_atms and atm not in grp1
                )
                possible_bonds = list(itertools.product(grp1, grp2))
                effective_bonds = (
                    bond for bond in possible_bonds if frozenset(bond) in gra[1].keys()
                )
                if len(tuple(set(effective_bonds).intersection(correct_bonds))) > 0:
                    grp = grp1 + grp2
                    if sorted(grp) not in [sorted(grpi) for grpi in grps]:
                        grps += (grp,)

    return grps
