"""Level 4 functions depending on other basic types (graph, geom)."""

from .. import amchi as amchi_
from .. import geom
from .. import graph as graph_
from ..extern import rdkit_
from .base import equivalent, has_stereo, standard_form


# # conversions
def graph(ich: str, stereo: bool = True):
    """Generate a molecular graph from an InChI string.

    :param ich: InChI string
    :param stereo: parameter to include stereochemistry information
    :rtype: automol molecular graph
    """
    gra = amchi_.graph(ich, stereo=stereo)
    if stereo:
        gra = graph_.with_explicit_stereo_hydrogens(gra)
    return gra


def geometry(ich, check=True):
    """Generate a molecular geometry from an InChI string.

    :param ich: InChI string
    :type ich: str
    :param check: check stereo and connectivity?
    :type check: bool
    :rtype: automol molecular geometry data structure
    """
    geo = amchi_.geometry(ich, check=check)
    return geo


def zmatrix(ich, check=True):
    """Generate a z-matrix from an InChI string.

    :param ich: InChI string
    :type ich: str
    :param check: check stereo and connectivity?
    :type check: bool
    :rtype: automol z-matrix data structure
    """
    zma = amchi_.zmatrix(ich, check=check)
    return zma


def amchi(ich, stereo=True):
    """Convert an InChI to an AMChI string.

    Only for good InChIs, where this can be done validly.

    :param ich: InChI string
    :type ich: str
    :returns: AMChI string
    :rtype: str
    """
    gra = graph(ich, stereo=stereo)

    if stereo:
        assert not is_bad(ich, gra=gra), (
            "Don't use this function with bad InChIs if stereo=True. "
            f"Bad InChI: {ich}"
        )

    ach = graph_.amchi(gra)
    return ach


def rdkit_molecule(ich, stereo=True):
    """Convert a InChI string to an RDKit molecule.

    This is mainly useful for quick visualization with IPython, which can
    be done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_molecule(ich))

    :param ich: InChI string
    :type ich: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :returns: the RDKit molecule
    """
    rdkit_.turn_3d_visualization_off()
    gra = graph(ich, stereo=stereo)
    return graph_.rdkit_molecule(gra, stereo=stereo)


def rdkit_reaction(richs, pichs, stereo=True, res_stereo=False):
    """Convert reactant and product graphs to an RDKit reaction object.

    This is mainly useful for quick visualization with IPython, which can be
    done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_reaction(pgras, rgras))

        :param richs: InChI strings for the reactants
        :param pichs: InChI strings for the products
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :returns: the RDKit reaction
    """
    rdkit_.turn_3d_visualization_off()
    rgras = [graph(s, stereo=stereo) for s in richs]
    pgras = [graph(s, stereo=stereo) for s in pichs]
    return graph_.rdkit_reaction(rgras, pgras, stereo=stereo, res_stereo=res_stereo)


def display(ich, stereo=True):
    """Display graph to IPython using the RDKit visualizer.

    :param ich: InChI string
    :type ich: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    """
    rdkit_.turn_3d_visualization_off()
    gra = graph(ich, stereo=stereo)
    graph_.display(gra, stereo=stereo)


def display_reaction(richs, pichs, stereo=True):
    """Display reaction to IPython using the RDKit visualizer.

    :param richs: InChI strings for the reactants
    :param pichs: InChI strings for the products
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    """
    rdkit_.turn_3d_visualization_off()
    rgras = [graph(s, stereo=stereo) for s in richs]
    pgras = [graph(s, stereo=stereo) for s in pichs]
    graph_.display_reaction(rgras, pgras, stereo=stereo)


# # derived properties
def is_complete(ich):
    """Determine if the InChI string is complete
    (has all stereo-centers assigned).

    :param ich: InChI string
    :type ich: str
    :rtype: bool
    """
    gra0 = graph(ich, stereo=False)

    needs_stereo = False
    if bool(graph_.unassigned_stereocenter_keys(gra0)):
        needs_stereo = not has_stereo(ich)

    return equivalent(ich, standard_form(ich)) and not needs_stereo


def is_bad(ich, gra=None):
    """Determine if the InChI string is bad, i.e. one of the InChI failure
    cases (resonance bond stereo, vinyl bond stereo, etc.

    :param ich: InChI string
    :type ich: str
    :param gra: A graph version of the InChI, to avoid recalculating
    :type gra: automol graph
    :returns: True if it is a bad InChI
    :rtype: bool
    """
    gra = graph(ich) if gra is None else gra
    ret = graph_.base.inchi_is_bad(gra, ich)
    return ret


# # derived transformations
def add_stereo(ich):
    """Add stereochemistry to an InChI string converting to/from geometry.

    :param ich: InChI string
    :type ich: str
    :rtype: str
    """
    geo = geometry(ich)
    ich = geom.inchi(geo, stereo=True)
    return ich


def expand_stereo(ich):
    """Obtain all possible stereoisomers of an InChI string.

    :param ich: InChI string
    :type ich: str
    :rtype: list[str]
    """
    gra = graph(ich, stereo=False)
    sgrs = graph_.expand_stereo(gra)
    ste_ichs = [graph_.inchi(sgr, stereo=True) for sgr in sgrs]
    return ste_ichs
