""" Level 4 functions depending on other basic types (geom, graph)
"""
import automol.graph
import automol.geom
from automol.extern import rdkit_
from automol.extern import pybel_
from automol.amchi.base import isotope_layers
from automol.amchi.base import symbols
from automol.amchi.base import bonds
from automol.amchi.base import hydrogen_valences
from automol.amchi.base import atom_stereo_parities
from automol.amchi.base import bond_stereo_parities
from automol.amchi.base import is_inverted_enantiomer
from automol.amchi.base import has_stereo
from automol.amchi.base import split


# # conversions
def smiles(chi, res_stereo=True):
    """ Convert a ChI string into a SMILES string.

        :param chi: ChI string
        :type chi: str
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :returns: the SMILES string
        :rtype: str
    """
    chis = split(chi)
    smis = [_connected_smiles(c, res_stereo=res_stereo) for c in chis]
    smi = '.'.join(smis)
    return smi


def _connected_smiles(chi, res_stereo=True):
    """ Convert a single-component ChI string into a SMILES string.

        :param chi: ChI string
        :type chi: str
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :returns: the SMILES string
        :rtype: str
    """
    gra = _connected_graph(chi, stereo=True, can=False)
    smi = automol.graph.smiles(gra, stereo=True, res_stereo=res_stereo)
    return smi


def graph(chi, stereo=True, can=False):
    """ Generate a molecular graph from a ChI string.

        :param chi: ChI string
        :type chi: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph
    """
    chis = split(chi)
    gras = [_connected_graph(c, stereo=stereo, can=can) for c in chis]
    gra = automol.graph.union_from_sequence(gras, shift_keys=True)
    return gra


def _connected_graph(chi, stereo=True, can=False):
    """ Generate a connected molecular graph from a single-component ChI string.

        :param chi: ChI string
        :type chi: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph
    """
    symb_dct = symbols(chi)
    bnd_keys = bonds(chi)
    atm_imp_hyd_vlc_dct = hydrogen_valences(chi)

    if isotope_layers(chi):
        raise NotImplementedError("Isotopic graph conversion not implemented")

    if stereo:
        bnd_ste_par_dct = bond_stereo_parities(chi)
        atm_ste_par_dct = atom_stereo_parities(chi)
    else:
        atm_ste_par_dct = None
        bnd_ste_par_dct = None

    is_inv = is_inverted_enantiomer(chi)

    gra = automol.graph.from_data(
        atm_symb_dct=symb_dct,
        bnd_keys=bnd_keys,
        atm_imp_hyd_vlc_dct=atm_imp_hyd_vlc_dct,
        atm_ste_par_dct=atm_ste_par_dct,
        bnd_ste_par_dct=bnd_ste_par_dct,
    )

    if is_inv is True:
        gra = automol.graph.reflect_local_stereo(gra)
        gra = automol.graph.from_local_stereo(gra)
    elif has_stereo(chi) and not can:
        gra = automol.graph.from_local_stereo(gra)

    return gra


# def geometry(ich, check=True):
#     """ Generate a molecular geometry from a ChI string
#     """
def _connected_geometry(chi, check=True):
    """ Generate a connected molecular geometry from a single-component ChI
        string.

        :param chi: ChI string
        :type chi: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular geometry
    """
    # Convert graph to local stereo to avoid multiple recanonicalizations
    gra = _connected_graph(chi, stereo=True)
    gra = automol.graph.explicit(gra)
    gra = automol.graph.to_local_stereo(gra)

    smi = _connected_smiles(chi, res_stereo=False)
    has_ste = has_stereo(chi)

    def _gen1():
        rdm = rdkit_.from_smiles(smi)
        geo, = rdkit_.to_conformers(rdm, nconfs=1)
        return geo

    def _gen2():
        pbm = pybel_.from_smiles(smi)
        geo = pybel_.to_geometry(pbm)
        return geo

    # for gen_ in (_gen1, _gen1, _gen1, _gen2):
    #     success = False
    #     try:

    geo = _gen1()

    # If the ChI had stereo, enforce correct stereo on the geometry. For
    # resonance bonds in particular, the stereo parity of the geometry is
    # arbitrary and needs to be corrected.
    if has_stereo(gra):
        gra_ = automol.geom.graph(geo)
        geo_idx_dct = automol.graph.isomorphism(gra, gra_, stereo=False)
        geo = automol.graph.stereo_corrected_geometry(
            gra, geo, geo_idx_dct=geo_idx_dct, loc=True)

    # Check to see if the geometry matches

    print(automol.graph.string(gra))
    print(automol.geom.string(geo))
    print(has_ste)
    print(geo_idx_dct)


if __name__ == '__main__':
    CHI = 'AMChI=1/C5H6FO/c6-4-2-1-3-5-7/h1-5,7H/b2-1-,3-1-,4-2-,5-3+'
    CHI = 'AMChI=1/C5H6FO/c6-4-2-1-3-5-7/h1-5,7H/b2-1+,3-1+,4-2+,5-3+'
    _connected_geometry(CHI, check=True)
