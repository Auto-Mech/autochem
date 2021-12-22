""" Level 4 functions depending on other basic types (geom, graph)
"""
import automol.graph.base
from automol.amchi.base import symbols
from automol.amchi.base import bonds
from automol.amchi.base import hydrogen_valences


# # conversions
def connected_graph(chi, stereo=True):
    """ Generate a molecular graph from a ChI string.

        :param chi: ChI string
        :type chi: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph
    """
    if stereo:
        raise NotImplementedError("this is in progress...")

    symb_dct = symbols(chi)
    bnds = bonds(chi)
    nhyd_dct = hydrogen_valences(chi)

    gra = automol.graph.base.from_data(
        atm_symb_dct=symb_dct,
        bnd_keys=bnds,
        atm_imp_hyd_vlc_dct=nhyd_dct
    )
    return gra


if __name__ == '__main__':
    GRA = ({0: ('C', 1, None), 1: ('C', 1, True), 2: ('C', 1, True),
            3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
            6: ('F', 0, None), 7: ('F', 0, None)},
           {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
            frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
            frozenset({2, 7}): (1, None)})
    print(automol.graph.string(GRA))
    CHI = automol.graph.amchi(GRA)
    print(CHI)
    GRA = connected_graph(CHI, stereo=False)
    print(automol.graph.string(GRA))
