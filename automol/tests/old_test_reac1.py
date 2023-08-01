""" test automol.reac
"""

import automol


def test__canonical_enantiomer():
    """ test reac.canonical_enantiomer
    """
    rct_smis = ['CC(OO)C(O[O])C(OO)C']
    prd_smis = ['CC(OO)C(OO)C(OO)[CH2]']

    rxn = automol.reac.with_structures_from_smiles(rct_smis, prd_smis)[0][0]

    # 2A. Full expansion -- includes non-canonical enantiomer reactions
    print("Full reaction expansion:")
    for srxn in automol.reac.expand_stereo(rxn, enant=True):
        rct_chis, prd_chis = automol.reac.chi(srxn)
        print(' +\n'.join(rct_chis) + " =>\n" + ' +\n'.join(prd_chis))

        # These functions operate directly on the reaction object:
        is_can = automol.reac.is_canonical_enantiomer(srxn)
        print(f"Canonical? {is_can}")
        # Convert it to a canonical enantiomer reaction like this
        srxn = automol.reac.canonical_enantiomer(srxn)
        assert automol.reac.is_canonical_enantiomer(srxn)

        # These are the equivalent functions for ChIs
        is_can = automol.chi.is_canonical_enantiomer_reaction(rct_chis,
                                                              prd_chis)
        print(f"Canonical? {is_can}")
        # Convert it to a canonical enantiomer reaction like this
        rct_chis, prd_chis = automol.chi.canonical_enantiomer_reaction(
            rct_chis, prd_chis)
        assert automol.chi.is_canonical_enantiomer_reaction(rct_chis, prd_chis)
        print()

    # 2B. Restricted expansion -- includes only canonical enantiomers
    print("Restricted reaction expansion:")
    for srxn in automol.reac.expand_stereo(rxn, enant=False):
        rct_chis, prd_chis = automol.reac.chi(srxn)
        print(' +\n'.join(rct_chis) + " =>\n" + ' +\n'.join(prd_chis))

        # Check canonicity for a reaction object
        assert automol.reac.is_canonical_enantiomer(srxn)

        # Check canonicity for reaction ChIs
        assert automol.chi.is_canonical_enantiomer_reaction(rct_chis, prd_chis)
        print()


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings("error")

    # test__add_stereo_from_unordered_geometries()
    # test__stereo()
    # test__canonical_enantiomer()
