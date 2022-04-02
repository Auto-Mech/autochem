""" functions operating on reaction ChI tuples
"""
from automol.amchi.base._core import sorted_
from automol.amchi.base._core import standard_form
from automol.amchi.base._core import is_chiral
from automol.amchi.base._core import reflect


def filter_enantiomer_reactions(rxn_chis_lst):
    """ Filter out mirror images from a list of reaction ChIs

        Redundant reactions are identified by inverting the enantiomers on both
        sides of the reaction (if reactants and/or products are enationmers)
        and removing them from the list.

        The list is sorted first, so that the same enantiomers will be chosen
        regardless of the order of their appearance.
    """
    rxn_chis_lst = sort_reactions(rxn_chis_lst)

    uniq_rxn_chis_lst = []
    seen_rxn_chis_lst = []

    for rxn_chis in rxn_chis_lst:
        rct_chis, prd_chis = rxn_chis

        # Standardize for comparison to seen ChIs
        rxn_chis = (tuple(map(standard_form, rct_chis)),
                    tuple(map(standard_form, prd_chis)))

        # Only proceed if we haven't seen the reaction before. A reaction in
        # the list of seen reactions will be ignored.
        if rxn_chis not in seen_rxn_chis_lst:
            uniq_rxn_chis_lst.append(rxn_chis)
            seen_rxn_chis_lst.append(rxn_chis)
            seen_rxn_chis_lst.append(tuple(reversed(rxn_chis)))

            # If the reaction is chiral, also add its mirror image to the list.
            if (any(map(is_chiral, rct_chis)) or
                    any(map(is_chiral, prd_chis))):
                refl_rxn_chis = (tuple(map(reflect, rct_chis)),
                                 tuple(map(reflect, prd_chis)))
                seen_rxn_chis_lst.append(refl_rxn_chis)
                seen_rxn_chis_lst.append(tuple(reversed(refl_rxn_chis)))

    uniq_rxn_chis_lst = tuple(
        (tuple(p), tuple(r)) for p, r in uniq_rxn_chis_lst)
    return uniq_rxn_chis_lst


def sort_reactions(rxn_chis_lst):
    """ Sort reactions such that enantiomeric versions always appear in the
        same order.
    """

    def _sort_key(rxn_chis):
        # Sort the reactants and products
        chis1, chis2 = map(sorted_, rxn_chis)

        # Reverse the reaction based on the autofile criterion
        if (len(chis1), chis1) < (len(chis2), chis2):
            chis1, chis2 = chis2, chis1

        # Return the fully sorted reaction
        return (chis1, chis2)

    rxn_chis_lst = sorted(rxn_chis_lst, key=_sort_key)
    return rxn_chis_lst