""" reaction formulae
"""
import itertools
from automol.formula._formula import add_hydrogen as _add_hydrogen


def argsort_hydrogen_abstraction(rct_fmls, prd_fmls):
    """ returns reactants and products sorted as RH + Q => R + QH
    """
    rxn_idxs = None
    if len(rct_fmls) == len(prd_fmls) == 2:
        for idxs1, idxs2 in itertools.product(
                itertools.permutations(range(2)), repeat=2):
            fmls1 = list(map(rct_fmls.__getitem__, idxs1))
            fmls2 = list(map(prd_fmls.__getitem__, idxs2))
            if (fmls1[0] == _add_hydrogen(fmls2[1]) and
                    _add_hydrogen(fmls1[1]) == fmls2[0]):
                rxn_idxs = (idxs1, idxs2)
                break
    return rxn_idxs
