"""reaction formulae."""
import itertools

from _collections_abc import Sequence

from ._form import Formula, add_element, join_sequence


def is_valid_reaction(rct_fmls: Sequence[Formula], prd_fmls: Sequence[Formula]) -> bool:
    """Use the formula to see if a reaction preserves stoichiometry.

    :param rct_fmls: stoichiometries of the reactants
    :param prd_fmls: stoichiometries of the products
    :return: True if the reaction is valid, False if not
    """
    return join_sequence(rct_fmls) == join_sequence(prd_fmls)


def argsort_hydrogen_abstraction(
    rct_fmls: Sequence[Formula], prd_fmls: Sequence[Formula]
) -> tuple[tuple[int, int], tuple[int, int]] | None:
    """Generate the indices which allows the reactants and products of
    a hydrogen abstraction reaction can be sorted as RH + Q => R + QH.

    :param rct_fmls: stoichiometries of the reactants
    :param prd_fmls: stoichiometries of the products
    :return: Indices to sort reaction to R + QH
    """
    rxn_idxs = None
    if len(rct_fmls) == len(prd_fmls) == 2:
        for idxs1, idxs2 in itertools.product(
            itertools.permutations(range(2)), repeat=2
        ):
            fmls1 = list(map(rct_fmls.__getitem__, idxs1))
            fmls2 = list(map(prd_fmls.__getitem__, idxs2))
            if (
                fmls1[0] == add_element(fmls2[1], "H")
                and add_element(fmls1[1], "H") == fmls2[0]
            ):
                rxn_idxs = (idxs1, idxs2)
                break

    return rxn_idxs
