""" functions operating on reaction ChI tuples
"""
from automol.amchi.base._core import split
from automol.amchi.base._core import reflect
from automol.amchi.base._core import sorted_
from automol.amchi.base._core import argsort
from automol.amchi.base._core import standard_form
from automol.amchi.base._core import is_chiral
from automol.amchi.base._core import is_canonical_enantiomer_list


def is_canonical_enantiomer_reaction(rct_chi, prd_chi):
    """ Does this reaction have a canonical combination of enantiomers?

        :param rct_chi: A multi-component ChI or list of ChIs for the reactants
        :type rct_chis: str or list[str]
        :param prd_chi: A multi-component ChI or list of ChIs for the products
        :type prd_chis: str or list[str]

        :returns: Whether or not the reaction is canonical
        :rtype: bool
    """
    rct_chis = split(rct_chi) if isinstance(rct_chi, str) else rct_chi
    prd_chis = split(prd_chi) if isinstance(prd_chi, str) else prd_chi

    # Switch to the canonical reaction direction
    if not is_canonical_reaction_direction(rct_chis, prd_chis):
        rct_chis, prd_chis = prd_chis, rct_chis

    chis = sorted_(rct_chis) + sorted_(prd_chis)
    can = is_canonical_enantiomer_list(chis)
    return can


def is_canonical_reaction_direction(rct_chis, prd_chis):
    """ Is this the canonical reaction direction, or should it be reversed?

        :param rct_chis: A list of ChIs for the reactants
        :type rct_chis: list[str]
        :param prd_chis: A list of ChIs for the products
        :type prd_chis: list[str]
        :returns: Whether or not the reaction is canonical
        :rtype: bool
    """
    nrcts = len(rct_chis)
    nprds = len(prd_chis)

    idxs = argsort(rct_chis + prd_chis)
    rct_idxs = sorted(idxs[:nrcts])
    prd_idxs = sorted(idxs[nrcts:])

    rct_rep = (nrcts, rct_idxs)
    prd_rep = (nprds, prd_idxs)
    return rct_rep < prd_rep


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
            if any(map(is_chiral, rct_chis + prd_chis)):
                refl_rxn_chis = (tuple(map(reflect, rct_chis)),
                                 tuple(map(reflect, prd_chis)))
                seen_rxn_chis_lst.append(refl_rxn_chis)
                seen_rxn_chis_lst.append(tuple(reversed(refl_rxn_chis)))

    uniq_rxn_chis_lst = tuple(uniq_rxn_chis_lst)
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
