""" functions operating on reaction InChI tuples
"""
from automol.inchi.base._core import sorted_
from automol.inchi.base._core import standard_form
from automol.inchi.base._core import is_enantiomer
from automol.inchi.base._core import reflect


def filter_enantiomer_reactions(rxn_ichs_lst):
    """ Filter out mirror images from a list of reaction InChIs

        Redundant reactions are identified by inverting the enantiomers on both
        sides of the reaction (if reactants and/or products are enationmers)
        and removing them from the list.

        The list is sorted first, so that the same enantiomers will be chosen
        regardless of the order of their appearance.
    """
    rxn_ichs_lst = sort_reactions(rxn_ichs_lst)

    uniq_rxn_ichs_lst = []
    seen_rxn_ichs_lst = []

    for rxn_ichs in rxn_ichs_lst:
        rct_ichs, prd_ichs = rxn_ichs

        # Standardize for comparison to seen InChIs
        rxn_ichs = (tuple(map(standard_form, rct_ichs)),
                    tuple(map(standard_form, prd_ichs)))

        # Only proceed if we haven't seen the reaction before. A reaction in
        # the list of seen reactions will be ignored.
        if rxn_ichs not in seen_rxn_ichs_lst:
            uniq_rxn_ichs_lst.append(rxn_ichs)
            seen_rxn_ichs_lst.append(rxn_ichs)
            seen_rxn_ichs_lst.append(tuple(reversed(rxn_ichs)))

            # If the reaction is chiral, also add its mirror image to the list.
            if (any(map(is_enantiomer, rct_ichs)) or
                    any(map(is_enantiomer, prd_ichs))):
                refl_rxn_ichs = (tuple(map(reflect, rct_ichs)),
                                 tuple(map(reflect, prd_ichs)))
                seen_rxn_ichs_lst.append(refl_rxn_ichs)
                seen_rxn_ichs_lst.append(tuple(reversed(refl_rxn_ichs)))

    uniq_rxn_ichs_lst = tuple(
        (tuple(p), tuple(r)) for p, r in uniq_rxn_ichs_lst)
    return uniq_rxn_ichs_lst


def sort_reactions(rxn_ichs_lst):
    """ Sort reactions such that enantiomeric versions always appear in the
        same order.
    """

    def _sort_key(rxn_ichs):
        # Sort the reactants and products
        ichs1, ichs2 = map(sorted_, rxn_ichs)

        # Reverse the reaction based on the autofile criterion
        if (len(ichs1), ichs1) < (len(ichs2), ichs2):
            ichs1, ichs2 = ichs2, ichs1

        # Return the fully sorted reaction
        return (ichs1, ichs2)

    rxn_ichs_lst = sorted(rxn_ichs_lst, key=_sort_key)
    return rxn_ichs_lst
