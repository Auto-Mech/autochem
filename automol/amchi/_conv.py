""" Level 4 functions depending on other basic types (geom, graph)
"""
import automol.graph.base
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
    gra = automol.graph.base.union_from_sequence(gras, shift_keys=True)
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

    gra = automol.graph.base.from_data(
        atm_symb_dct=symb_dct,
        bnd_keys=bnd_keys,
        atm_imp_hyd_vlc_dct=atm_imp_hyd_vlc_dct,
        atm_ste_par_dct=atm_ste_par_dct,
        bnd_ste_par_dct=bnd_ste_par_dct,
    )

    if is_inv is True:
        gra = automol.graph.base.reflect_local_stereo(gra)
        gra = automol.graph.base.from_local_stereo(gra)
    elif has_stereo(chi) and not can:
        gra = automol.graph.base.from_local_stereo(gra)

    return gra


if __name__ == '__main__':
    GRA = ({0: ('C', 1, None), 1: ('C', 1, True), 2: ('C', 1, True),
            3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
            6: ('F', 0, None), 7: ('F', 0, None)},
           {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
            frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
            frozenset({2, 7}): (1, None)})
    GRA = ({0: ('C', 1, False), 1: ('C', 1, True), 2: ('C', 1, False),
            3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
            6: ('F', 0, None), 7: ('F', 0, None)},
           {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
            frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
            frozenset({2, 7}): (1, None)})
    print(automol.graph.string(GRA))
    CHI = automol.graph.amchi(GRA)
    print(CHI)
    GRA = _connected_graph(CHI)
    print(automol.graph.string(GRA))
