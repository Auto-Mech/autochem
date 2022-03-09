""" Level 4 functions depending on other basic types (geom, graph)
"""
import automol.graph.base
from automol.rsmiles.base import parse_properties


# # conversions
def connected_graph(smi, stereo=True):
    """ Generate a molecular graph from a ChI string.

        :param smi: SMILES string
        :type smi: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph
    """
    symb_dct, bnd_ord_dct, atm_par_dct, bnd_par_dct = parse_properties(smi)
    bnd_keys = bnd_ord_dct.keys()

    if not stereo:
        atm_par_dct = None
        bnd_par_dct = None

    gra = automol.graph.base.from_data(
        atm_symb_dct=symb_dct,
        bnd_keys=bnd_keys,
        atm_ste_par_dct=atm_par_dct,
        bnd_ste_par_dct=bnd_par_dct,
    )

    if automol.graph.base.has_stereo(gra):
        # The parser marks all bonds with directional bonds on either side as
        # having stereo, because it has no way to distinguish between them.  In
        # lieu of a more rigorous check, remove stereo from all non-sp2 bonds.
        # If this is an issue, we could create a more rigorous check to see if
        # a bond is stereogenic.
        ste_bnd_keys = automol.graph.base.bond_stereo_keys(gra)
        sp2_bnd_keys = automol.graph.base.sp2_bond_keys(gra)
        bnd_keys = ste_bnd_keys - sp2_bnd_keys
        gra = automol.graph.base.remove_bond_stereo_parities(gra, bnd_keys)

        # Convert from local to canonical stereo
        gra = automol.graph.base.from_local_stereo(gra)

    return gra
