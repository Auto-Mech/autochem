""" Level 4 functions depending on amchi, inchi, graph, geom
"""
import automol.amchi
import automol.inchi


# # conversions
# # derived properties
def is_complete(chi):
    """ Determine if the ChI string is complete
        (has all stereo-centers assigned).

        Currently only checks species that does not have any
        resonance structures.

        :param chi: ChI string
        :type chi: str
        :rtype: bool
    """
    pfx = automol.amchi.base.prefix(chi)
    if pfx == 'AMChI':
        ret = True
    elif pfx == 'InChI':
        ret = automol.inchi.is_compoete(chi)
    else:
        raise ValueError(f"ChI string '{chi}' has unknown prefix '{pfx}'.")
    return ret


# # derived transformations
def add_stereo(chi):
    """ Add stereochemistry to a ChI string converting to/from geometry.

        :param chi: ChI string
        :type chi: str
        :rtype: str
    """
    pfx = automol.amchi.base.prefix(chi)
    if pfx == 'AMChI':
        ret = automol.amchi.add_stereo(chi)
    elif pfx == 'InChI':
        ret = automol.inchi.add_stereo(chi)
    else:
        raise ValueError(f"ChI string '{chi}' has unknown prefix '{pfx}'.")
    return ret


def expand_stereo(chi):
    """ Obtain all possible stereoisomers of a ChI string.

        :param chi: ChI string
        :type chi: str
        :rtype: list[str]
    """
    gra = automol.amchi.graph(chi, stereo=False)
    sgrs = automol.graph.stereomers(gra)
    ste_chis = []
    for sgr in sgrs:
        # If the graph has resonance bond stereo, use an AMChI string
        if automol.graph.has_resonance_bond_stereo(sgr):
            ste_chi = automol.graph.amchi(sgr, stereo=True)
        else:
            ste_chi = automol.graph.inchi(sgr, stereo=True)

        ste_chis.append(ste_chi)

    return ste_chis
