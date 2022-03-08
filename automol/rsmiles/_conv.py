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
    symb_dct, atm_imp_hyd_vlc_dct, bnd_ord_dct = parse_properties(smi)
    bnd_keys = bnd_ord_dct.keys()

    if stereo:
        raise NotImplementedError("Not yet implemented.")
    else:
        atm_ste_par_dct = None
        bnd_ste_par_dct = None

    gra = automol.graph.base.from_data(
        atm_symb_dct=symb_dct,
        bnd_keys=bnd_keys,
        atm_imp_hyd_vlc_dct=atm_imp_hyd_vlc_dct,
        atm_ste_par_dct=atm_ste_par_dct,
        bnd_ste_par_dct=bnd_ste_par_dct,
    )

    return gra
