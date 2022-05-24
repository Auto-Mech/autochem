""" Level 4 functions depending on other basic types (geom, graph)
"""
import functools
import automol.graph
import automol.geom
from automol import error
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
from automol.amchi.base import with_inchi_prefix
from automol.amchi.base import equivalent
from automol.amchi.base import standard_form


# # conversions
def amchi_key(chi):
    """ Generate a ChIKey from a ChI string.

        :param chi: ChI string
        :type chi: str
        :rtype: str
    """
    ich = with_inchi_prefix(chi)
    return rdkit_.inchi_to_inchi_key(ich)


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
    gra = _connected_graph(chi, stereo=True, local_stereo=True)
    smi = automol.graph.smiles(gra, stereo=True, local_stereo=True,
                               res_stereo=res_stereo)
    return smi


def graph(chi, stereo=True, local_stereo=False):
    """ Generate a molecular graph from a ChI string.

        :param chi: ChI string
        :type chi: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :param local_stereo: assign local stereo parities?
        :type local_stereo: bool
        :rtype: automol molecular graph
    """
    chis = split(chi)
    gras = [_connected_graph(c, stereo=stereo, local_stereo=local_stereo)
            for c in chis]
    gra = automol.graph.union_from_sequence(gras, shift_keys=True)
    return gra


def _connected_graph(chi, stereo=True, local_stereo=False):
    """ Generate a connected molecular graph from a single-component ChI string.

        :param chi: ChI string
        :type chi: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :param local_stereo: assign local stereo parities?
        :type local_stereo: bool
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

    if has_stereo(chi) and not local_stereo:
        gra = automol.graph.from_local_stereo(gra)

    return gra


def geometry(chi, check=True):
    """ Generate a molecular geometry from a ChI string.

        :param chi: ChI string
        :type chi: str
        :param check: check stereo and connectivity?
        :type check: bool
        :rtype: automol molecular geometry data structure
    """

    # rdkit fails for multi-component molecules, so we split it up and space
    # out the geometries
    chis = split(chi)
    geos = [_connected_geometry(chi, check=check) for chi in chis]
    geos = [automol.geom.translate(geo, [50. * idx, 0., 0.])
            for idx, geo in enumerate(geos)]
    geo = functools.reduce(automol.geom.join, geos)
    return geo


def _connected_geometry(chi, check=True):
    """ Generate a connected molecular geometry from a single-component ChI
        string.

        :param chi: ChI string
        :type chi: str
        :param check: check stereo and connectivity?
        :type check: bool
        :rtype: automol molecular geometry
    """
    # Convert graph to local stereo to avoid multiple recanonicalizations
    gra = _connected_graph(chi, stereo=True, local_stereo=True)
    gra = automol.graph.explicit(gra)

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

    def _gen3():
        if has_stereo(chi):
            raise ValueError

        gra = graph(chi, stereo=False)
        gra = automol.graph.explicit(gra)
        geo = automol.graph.embed.geometry(gra)
        return geo

    success = False
    for gen_ in (_gen1, _gen1, _gen1, _gen2, _gen3):
        try:
            geo = gen_()
        except (RuntimeError, TypeError, ValueError):
            continue

        # If the ChI has stereo, enforce correct stereo on the geometry.
        if check:
            # There is stereo.
            # First, check connectivity.
            gra_ = automol.geom.graph(geo)
            geo_idx_dct = automol.graph.isomorphism(
                    gra, gra_, stereo=False)

            if geo_idx_dct is None:
                continue

            geo = automol.graph.linear_vinyl_corrected_geometry(
                gra, geo, geo_idx_dct=geo_idx_dct)

            if not has_ste:
                success = True
                break

            # Enforce correct stereo parities. This is necessary for
            # resonance bond stereo.
            geo = automol.graph.stereo_corrected_geometry(
                gra, geo, geo_idx_dct=geo_idx_dct, local_stereo=True)

            success = True
            break

    if not success:
        raise error.FailedGeometryGenerationError('Failed AMChI:', chi)

    return geo


def conformers(chi, nconfs=1, check=True, accept_fewer=False):
    """ Generate a connected molecular geometry from a single-component ChI
        string.

        :param chi: ChI string
        :type chi: str
        :param nconfs: number of conformer structures to generate
        :type nconfs: int
        :param check: check stereo and connectivity?
        :type check: bool
        :param accept_fewer: accept fewer than nconfs conformers?
        :type accept_fewer: bool
        :rtype: automol molecular geometry
    """
    # Convert graph to local stereo to avoid multiple recanonicalizations
    gra = _connected_graph(chi, stereo=True, local_stereo=True)
    gra = automol.graph.explicit(gra)

    smi = _connected_smiles(chi, res_stereo=False)
    has_ste = has_stereo(chi)

    try:
        rdm = rdkit_.from_smiles(smi)
        geos = rdkit_.to_conformers(rdm, nconfs=nconfs)
    except (RuntimeError, TypeError, ValueError) as err:
        raise error.FailedGeometryGenerationError(
            'Failed AMChI:', chi) from err

    ret_geos = []

    # If the ChI has stereo, enforce correct stereo on the geometry.
    if check:
        for geo in geos:
            if not has_ste:
                # If there wasn't stereo, only check connectivity
                gra_ = automol.geom.graph(geo, stereo=False)
                if automol.graph.isomorphism(gra, gra_, stereo=False):
                    ret_geos.append(geo)
            else:
                # There is stereo.
                # First, check connectivity.
                gra_ = automol.geom.graph(geo)
                geo_idx_dct = automol.graph.isomorphism(
                        gra, gra_, stereo=False)

                if geo_idx_dct is not None:
                    # Enforce correct stereo parities. This is necessary for
                    # resonance bond stereo.
                    geo = automol.graph.stereo_corrected_geometry(
                        gra, geo, geo_idx_dct=geo_idx_dct, local_stereo=True)

                    # Check if the assignment worked.
                    gra_ = automol.graph.set_stereo_from_geometry(gra_, geo)
                    if automol.graph.isomorphism(gra, gra_):
                        ret_geos.append(geo)

    if len(ret_geos) < nconfs and not accept_fewer:
        raise error.FailedGeometryGenerationError('Failed AMChI:', chi)

    return ret_geos


def zmatrix(chi, check=True):
    """ Generate a z-matrix from an InChI string.

        :param chi: InChI string
        :type chi: str
        :param check: check stereo and connectivity?
        :type check: bool
        :rtype: automol z-matrix data structure
    """

    geo = geometry(chi, check=check)
    zma = automol.geom.zmatrix(geo)
    return zma


# # derived properties
def is_complete(ich):
    """ Determine if the InChI string is complete
        (has all stereo-centers assigned).

        Currently only checks species that does not have any
        resonance structures.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """

    gra = graph(ich, stereo=False)
    ste_atm_keys = automol.graph.stereogenic_atom_keys(gra)
    ste_bnd_keys = automol.graph.stereogenic_bond_keys(gra)
    graph_has_stereo = bool(ste_atm_keys or ste_bnd_keys)

    _complete = equivalent(ich, standard_form(ich)) and not (
        has_stereo(ich) ^ graph_has_stereo)

    return _complete


# # derived transformations
def add_stereo(chi):
    """ Add stereochemistry to a ChI string converting to/from geometry.

        :param chi: ChI string
        :type chi: str
        :rtype: str
    """
    geo = geometry(chi)
    chi = automol.geom.amchi(geo, stereo=True)
    return chi


def expand_stereo(chi):
    """ Obtain all possible stereoisomers of a ChI string.

        :param chi: ChI string
        :type chi: str
        :rtype: list[str]
    """
    gra = graph(chi, stereo=False)
    sgrs = automol.graph.stereomers(gra)
    ste_chis = [automol.graph.amchi(sgr, stereo=True) for sgr in sgrs]
    return ste_chis


if __name__ == '__main__':
    CHI = 'AMChI=1/C5H6FO/c6-4-2-1-3-5-7/h1-5,7H/b2-1-,3-1-,4-2-,5-3+'
    CHI = 'AMChI=1/C5H6FO/c6-4-2-1-3-5-7/h1-5,7H/b2-1+,3-1+,4-2+,5-3+'
    # GEO = geometry(CHI, check=True)
    # print(automol.geom.string(GEO))

    for GEO in conformers(CHI, nconfs=5):
        print(automol.geom.string(GEO))
        print()
