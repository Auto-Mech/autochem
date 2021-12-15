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
    import numpy

    ACH = ('AMChI=1/C10H14ClFO/c1-8(9(5-12)10(13)6-11)7-3-2-4-7/'
           'h2-4,8-10,13H,5-6H2,1H3')
    GRA = connected_graph(ACH, stereo=False)

    NATMS = len(automol.graph.base.atom_keys(GRA))

    for _ in range(10):
        PMT = list(map(int, numpy.random.permutation(NATMS)))
        PMT_GRA = automol.graph.base.relabel(GRA, dict(enumerate(PMT)))
        PMT_ACH = automol.graph.base.amchi(PMT_GRA, stereo=False)
        print(automol.graph.base.string(PMT_GRA))
        print(PMT_ACH)
        assert PMT_ACH == ACH
