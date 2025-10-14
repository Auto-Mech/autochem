""" Level 4 functions depending on amchi, inchi, graph, geom
"""
from .. import amchi, geom, graph, inchi
from .base import reflect, split


# # conversions
# # derived properties
def is_complete(chi):
    """Determine if the ChI string is complete
    (has all stereo-centers assigned).

    Currently only checks species that does not have any
    resonance structures.

    :param chi: ChI string
    :type chi: str
    :rtype: bool
    """
    pfx = amchi.prefix(chi)
    if pfx == "AMChI":
        ret = amchi.is_complete(chi)
    elif pfx == "InChI":
        ret = inchi.is_complete(chi)
    else:
        raise ValueError(f"ChI string '{chi}' has unknown prefix '{pfx}'.")
    return ret


def is_bad(chi, gra=None):
    """Is this ChI string is a bad InChI?

    :param chi: ChI string
    :type chi: str
    :param gra: A graph version of the InChI, to avoid recalculating
    :type gra: automol graph
    :returns: True if it is a bad InChI
    :rtype: bool
    """
    pfx = amchi.prefix(chi)
    ret = False
    if pfx == "InChI":
        gra = inchi.graph(chi) if gra is None else gra
        ret = graph.inchi_is_bad(gra, chi)

    return ret


def chi_(chi):  # => chi() in __init__.py (underscore here avoids name clash)
    """Recalculate the ChI string, returning an AMChI only if the InChI is bad

    :param chi: ChI string
    :type chi: str
    :returns: an InChI if valid; otherwise an AMChI
    """
    gra = amchi.graph(chi)
    chi = graph.inchi(gra)
    if graph.inchi_is_bad(gra, chi):
        chi = graph.amchi(gra)
    return chi


def without_stereo(chi, reassess_amchi=True):
    """Remove stereo information from this ChI

    :param chi: ChI string
    :type chi: str
    """
    pfx = amchi.prefix(chi)
    if pfx == "AMChI":
        ret = amchi.without_stereo(chi)
    elif pfx == "InChI":
        ret = inchi.without_stereo(chi)
    else:
        raise ValueError(f"ChI string '{chi}' has unknown prefix '{pfx}'.")

    if reassess_amchi:
        ret = chi_(ret)

    return ret


# # derived transformations
def join(chis):
    """Join separate ChI strings into one multi-component ChI string.

    :param chis: sequence of ChI strings
    :type chis: tuple[str]
    :returns: the joined ChI string
    :rtype: str
    """
    pfxs = list(map(amchi.prefix, chis))
    if "AMChI" in pfxs:
        # If any of these are AMChIs, they must all be converted over before
        # joining them together.
        chis = [inchi.amchi(c) if amchi.prefix(c) == "InChI" else c for c in chis]
        ret = amchi.join(chis)
    elif set(pfxs) == {"InChI"}:
        ret = inchi.join(chis)
    else:
        raise ValueError(f"One of these has an unknown prefix: {chis}.")
    return ret


def add_stereo(chi):
    """Add stereochemistry to a ChI string converting to/from geometry.

    :param chi: ChI string
    :type chi: str
    :rtype: str
    """

    ste_slyrs = amchi.stereo_layers(chi)
    if not ste_slyrs:
        geo = amchi.geometry(chi)
        chi = geom.chi(geo, stereo=True)

    return chi


def expand_stereo(chi, enant=True):
    """Obtain all possible stereoisomers of a ChI string.

    :param chi: ChI string
    :type chi: str
    :param enant: Include all enantiomers, or only canonical ones?
    :type enant: bool
    :rtype: list[str]
    """
    gra = amchi.graph(chi, stereo=False)
    sgrs = graph.expand_stereo(gra, enant=enant, symeq=False)
    ste_chis = []
    for sgr in sgrs:
        ste_chi = graph.inchi(sgr, stereo=True)

        # Switch to AMChI if the InChI is bad
        if is_bad(ste_chi, gra=sgr):
            ste_chi = graph.amchi(sgr, stereo=True)

        ste_chis.append(ste_chi)

    return ste_chis


def is_canonical_enantiomer(chi):
    """Is this a canonical enantiomer? Also returns true for non-enantiomers.

    For multi-component InChIs, it checks that the first enantiomeric
    component (if any), in sorted order, is canonical.

    :param chi: ChI string
    :type chi: str
    :returns: False if the ChI is the inverted form of its canonical
        enantiomer; True in all other cases (including non-enantiomers).
    :rtype: bool
    """
    pfx = amchi.prefix(chi)
    if pfx == "InChI":
        chi = inchi.amchi(chi) if isinstance(chi, str) else tuple(map(inchi.amchi, chi))
    return amchi.is_canonical_enantiomer(chi)


def is_canonical_enantiomer_reaction(rct_chi, prd_chi):
    """Does this reaction have a canonical combination of enantiomers?

    :param rct_chi: A multi-component ChI or list of ChIs for the reactants
    :type rct_chi: str or list[str]
    :param prd_chi: A multi-component ChI or list of ChIs for the products
    :type prd_chi: str or list[str]

    :returns: Whether or not the reaction is canonical
    :rtype: bool
    """
    rct_chis = split(rct_chi) if isinstance(rct_chi, str) else rct_chi
    prd_chis = split(prd_chi) if isinstance(prd_chi, str) else prd_chi

    rct_chis = tuple(
        inchi.amchi(c) if amchi.prefix(c) == "InChI" else c for c in rct_chis
    )
    prd_chis = tuple(
        inchi.amchi(c) if amchi.prefix(c) == "InChI" else c for c in prd_chis
    )

    return amchi.is_canonical_enantiomer_reaction(rct_chis, prd_chis)


def canonical_enantiomer(chi):
    """Convert this ChI string to a canonical enantiomer, if it isn't already

    Works for multi-component InChIs or lists of InChIs

    :param chi: ChI string
    :type chi: str
    :returns: The reflected ChI string, if it was a non-canonical
        enantiomer; otherwise, it returns the original ChI string
    """
    if not is_canonical_enantiomer(chi):
        chi = reflect(chi) if isinstance(chi, str) else tuple(map(reflect, chi))
    return chi


def canonical_enantiomer_reaction(rct_chi, prd_chi):
    """Convert this reaction into a canonical combination of enantiomers

    :param rct_chi: A multi-component ChI or list of ChIs for the reactants
    :type rct_chi: str or list[str]
    :param prd_chi: A multi-component ChI or list of ChIs for the products
    :type prd_chi: str or list[str]

    :returns: The reflected reaction, if it was a non-canonical combination
        of enantiomers; otherwise, it returns the original reaction.
    """
    if not is_canonical_enantiomer_reaction(rct_chi, prd_chi):
        rct_chi = (
            reflect(rct_chi)
            if isinstance(rct_chi, str)
            else tuple(map(reflect, rct_chi))
        )
        prd_chi = (
            reflect(prd_chi)
            if isinstance(prd_chi, str)
            else tuple(map(reflect, prd_chi))
        )
    return rct_chi, prd_chi


# # InCHI string converters
def inchi_to_amchi(ich, convert=True, printlog=True):
    """For an InChI string with stereo convert to an AMChI
        Since stereo conversions are expensive, function checks
        if input string is an InChI, w/ stereo, for species
        with resonance stabilization. Then it attempts the
        conversion.
        If convert=True: converts bad inchis in any case

    Note that these results still depend on the InChI code
    algorithm not building an incorrect geometry, so this
    conversion is not guaranteed to work.
    """

    chi = ich
    if amchi.prefix(ich) == "InChI":
        if amchi.stereo_layers(ich) or convert:
            if is_bad(ich):
                geo = amchi.geometry(ich)
                chi = geom.chi(geo)

                # Print a warning since conversion not guaranteed
                if printlog:
                    if amchi.prefix(chi) == "AMChI":
                        print(
                            "Check accuracy of InChI -> AMChI conversion:"
                            f"{ich} -> {chi}"
                        )
    elif convert and amchi.prefix(chi) == "AMChI":
        gra = amchi.graph(chi)
        chi = graph.chi(gra)

    return chi
