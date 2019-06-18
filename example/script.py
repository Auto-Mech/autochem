""" sample script
"""
import os
import itertools
import pandas
import automol


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

    RCTS_ICH_LST = list(HEPTANE_TAB['reac_inchi'])
    PRDS_ICH_LST = list(HEPTANE_TAB['prod_inchi'])

    for IDX, (RCTS_ICH, PRDS_ICH) in enumerate(
            zip(RCTS_ICH_LST, PRDS_ICH_LST)):
        print(IDX)
        print(RCTS_ICH)
        print(PRDS_ICH)
        RCTS_CGR = automol.graph.explicit(automol.inchi.graph(RCTS_ICH))
        PRDS_CGR = automol.graph.explicit(automol.inchi.graph(PRDS_ICH))

        RET = classify(RCTS_CGR, PRDS_CGR)
        if RET is not None:
            TYP, RXN = RET
            print(TYP, RXN)
            assert automol.graph.backbone_isomorphic(
                automol.graph.reaction.react(RXN, RCTS_CGR), PRDS_CGR)

            for RCTS_SGR, PRDS_SGR in itertools.product(
                    automol.graph.stereomers(RCTS_CGR),
                    automol.graph.stereomers(PRDS_CGR)):
                print(automol.graph.reaction.is_stereo_compatible(
                    RXN, RCTS_SGR, PRDS_SGR))
        else:
            print('unclassified')
        print()
