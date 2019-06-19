""" sample script
"""
import os
import itertools
import pandas
import automol
import autofile


class ReactionType():
    """ reaction types """

    H_MIGRATION = 'HMIG'
    BETA_SCISSION = 'BSC'
    ADDITION = 'ADD'
    H_ABSTRACTION = 'HABS'


def classify(xgr1, xgr2):
    """ classify a reaction by type
    """
    ret = None

    rxn = automol.graph.reaction.hydrogen_migration(xgr1, xgr2)
    if rxn and ret is None:
        typ = ReactionType.H_MIGRATION
        ret = (typ, rxn)

    rxn = automol.graph.reaction.beta_scission(xgr1, xgr2)
    if rxn and ret is None:
        typ = ReactionType.BETA_SCISSION
        ret = (typ, rxn)

    rxn = automol.graph.reaction.addition(xgr1, xgr2)
    if rxn and ret is None:
        typ = ReactionType.ADDITION
        ret = (typ, rxn)

    rxn = automol.graph.reaction.hydrogen_abstraction(xgr1, xgr2)
    if rxn and ret is None:
        typ = ReactionType.H_ABSTRACTION
        ret = (typ, rxn)

    return ret


if __name__ == '__main__':
    PATH = os.path.dirname(os.path.realpath(__file__))
    HEPTANE_TAB = pandas.read_csv(os.path.join(PATH, 'data', 'heptane.csv'))

    ICH1_LST = list(HEPTANE_TAB['reac_inchi'])
    ICH2_LST = list(HEPTANE_TAB['prod_inchi'])
    MLT1_LST = list(HEPTANE_TAB['reac_mults'])
    MLT2_LST = list(HEPTANE_TAB['prod_mults'])

    for IDX, (ICH1, ICH2, MLT1, MLT2) in enumerate(
            zip(ICH1_LST, ICH2_LST, MLT1_LST, MLT2_LST)):
        ICHS1 = automol.inchi.split(ICH1)
        ICHS2 = automol.inchi.split(ICH2)
        MULTS1 = list(map(int, MLT1.split('_')))
        MULTS2 = list(map(int, MLT2.split('_')))
        ICHS_PAIR = (ICHS1, ICHS2)
        MULTS_PAIR = (MULTS1, MULTS2)
        CHARS_PAIR = ((0,) * len(ICHS1), (0,) * len(ICHS2))
        ICHS_PAIR, MULTS_PAIR, CHARS_PAIR = autofile.system.map_.sort_together(
            ICHS_PAIR, MULTS_PAIR, CHARS_PAIR)

        ICH1 = automol.inchi.join(ICHS_PAIR[0])
        ICH2 = automol.inchi.join(ICHS_PAIR[1])
        CGR1 = automol.graph.explicit(automol.inchi.graph(ICH1))
        CGR2 = automol.graph.explicit(automol.inchi.graph(ICH2))

        RET = classify(CGR1, CGR2)
        if RET is not None:
            TYP, RXN = RET
            print(IDX, TYP, RXN)
            for SGR1, SGR2 in itertools.product(
                    automol.graph.stereomers(CGR1),
                    automol.graph.stereomers(CGR2)):
                if automol.graph.reaction.is_stereo_compatible(
                        RXN, SGR1, SGR2):
                    SICH1 = automol.graph.inchi(SGR1)
                    SICH2 = automol.graph.inchi(SGR2)
                    SICHS1 = automol.inchi.split(SICH1)
                    SICHS2 = automol.inchi.split(SICH2)
                    SICHS_PAIR = (SICHS1, SICHS2)
                    print(SICHS_PAIR, CHARS_PAIR, MULTS_PAIR)
            print()
