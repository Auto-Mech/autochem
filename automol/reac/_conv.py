""" Conversion functions
"""
import IPython
import automol.graph
from automol.extern import rdkit_
from automol.reac._reac import reactant_graphs, product_graphs


# Conversion stuff
def amchi(rxn, stereo=True):
    """ Convert the reaction object to AMChIs

        :param rxn: the reaction object
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :returns: AMChI strings for the reactants and products
        :rtype: (tuple[str], tuple[str])
    """
    rct_chis = tuple(automol.graph.amchi(gra, stereo=stereo)
                     for gra in reactant_graphs(rxn))
    prd_chis = tuple(automol.graph.amchi(gra, stereo=stereo)
                     for gra in product_graphs(rxn))
    return (rct_chis, prd_chis)


def inchi(rxn, stereo=True):
    """ Convert the reaction object to InChIs

        :param rxn: the reaction object
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :returns: InChI strings for the reactants and products
        :rtype: (tuple[str], tuple[str])
    """
    rct_ichs = tuple(automol.graph.inchi(gra, stereo=stereo)
                     for gra in reactant_graphs(rxn))
    prd_ichs = tuple(automol.graph.inchi(gra, stereo=stereo)
                     for gra in product_graphs(rxn))
    return (rct_ichs, prd_ichs)


def chi(rxn, stereo=True):
    """ Convert the reaction object to ChIs

        :param rxn: the reaction object
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :returns: ChI strings for the reactants and products
        :rtype: (tuple[str], tuple[str])
    """
    rct_chis = tuple(automol.graph.chi(gra, stereo=stereo)
                     for gra in reactant_graphs(rxn))
    prd_chis = tuple(automol.graph.chi(gra, stereo=stereo)
                     for gra in product_graphs(rxn))
    return (rct_chis, prd_chis)


def smiles(rxn, stereo=True, res_stereo=True, exp_singles=False):
    """ Convert the reaction object to SMILESs

        :param rxn: the reaction object
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :param exp_singles: Use explicit '-' for single bonds?
        :type exp_singles: bool
        :returns: SMILES strings for the reactants and products
        :rtype: (tuple[str], tuple[str])
    """
    rct_smis = tuple(
        automol.graph.smiles(gra, stereo=stereo, res_stereo=res_stereo,
                             exp_singles=exp_singles)
        for gra in reactant_graphs(rxn))
    prd_smis = tuple(
        automol.graph.smiles(gra, stereo=stereo, res_stereo=res_stereo,
                             exp_singles=exp_singles)
        for gra in product_graphs(rxn))
    return (rct_smis, prd_smis)


def rdkit_reaction(rxn, stereo=True, res_stereo=False):
    """ Convert a reaction object to an RDKit reaction object.

    This is mainly useful for quick visualization with IPython, which can be
    done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_reaction(gra))

        :param rxn: the reaction object
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :param rxn: the reaction object
        :returns: the RDKit reaction
    """
    rdkit_.turn_3d_visualization_off()
    rct_smis, prd_smis = smiles(rxn, stereo=stereo, res_stereo=res_stereo,
                                exp_singles=True)
    rcts_smi = '.'.join(rct_smis)
    prds_smi = '.'.join(prd_smis)
    rxn_smi = f'{rcts_smi}>>{prds_smi}'
    return rdkit_.from_smarts(rxn_smi)


def display(rxn):
    """ Display reaction object to IPython using the RDKit visualizer

        :param rxn: the reaction object
        :returns: None
    """
    rdkit_.turn_3d_visualization_off()
    return IPython.display.display(rdkit_reaction(rxn))


# Get a reaction object from various identifiers
def with_structures_from_chi(rct_chis, prd_chis, zmat=False, stereo=False,
                             ts_stereo=None, ts_enant=None):
    """ Get reaction objects with geometry/z-matrix structures from ChIs

    :param rct_chis: The reactant ChI (InChI or AMChI) strings
    :type rct_chis: list[str]
    :param prd_chis: The product ChI (InChI or AMChI) strings
    :type prd_chis: list[str]
    :param zmat: Return z-matrix structures instead of cartesian geometries?
    :type zmat: bool
    :param stereo: Include stereoassignments? If `True`, then `ts_stereo` will
        default to `True` and `ts_enant` will default to `False`.
    :type stereo: bool
    :param ts_stereo: Include fleeting TS stereoisomers as distinct TSs?
    :type ts_stereo: bool or NoneType
    :param ts_enant: Include fleeting TS enantiomers as distinct TSs?
    :type ts_enant: bool or NoneType
    :returns: A series of tuples containing, in order, the reaction object, the
        TS structure, the (sorted) reactant structures, and the (sorted)
        product structures. The latter three are returned as either automol
        geometry objects or, if `zmatrix` was set to `True, automol z-matrix
        objects.
    """

    rct_geos = list(map(automol.chi.geometry, rct_chis))
    prd_geos = list(map(automol.chi.geometry, prd_chis))

    return with_structures_from_geometry(
        rct_geos, prd_geos, zmat=zmat, stereo=stereo, ts_stereo=ts_stereo,
        ts_enant=ts_enant)


def with_structures_from_smiles(rct_smis, prd_smis, zmat=False, stereo=False,
                                ts_stereo=None, ts_enant=None):
    """ Get reaction objects with geometry/z-matrix structures from SMILES

    :param rct_smis: The reactant SMILES strings
    :type rct_smis: list[str]
    :param prd_smis: The product SMILES strings
    :type prd_smis: list[str]
    :param zmat: Return z-matrix structures instead of cartesian geometries?
    :type zmat: bool
    :param stereo: Include stereoassignments? If `True`, then `ts_stereo` will
        default to `True` and `ts_enant` will default to `False`.
    :type stereo: bool
    :param ts_stereo: Include fleeting TS stereoisomers as distinct TSs?
    :type ts_stereo: bool or NoneType
    :param ts_enant: Include fleeting TS enantiomers as distinct TSs?
    :type ts_enant: bool or NoneType
    :returns: A series of tuples containing, in order, the reaction object, the
        TS structure, the (sorted) reactant structures, and the (sorted)
        product structures. The latter three are returned as either automol
        geometry objects or, if `zmatrix` was set to `True, automol z-matrix
        objects.
    """

    # rct_geos = list(map(automol.smiles.geometry, rct_smis))
    # prd_geos = list(map(automol.smiles.geometry, prd_smis))

    rct_chis = list(map(automol.smiles.chi, rct_smis))
    prd_chis = list(map(automol.smiles.chi, prd_smis))
    rct_geos = list(map(automol.chi.geometry, rct_chis))
    prd_geos = list(map(automol.chi.geometry, prd_chis))

    return with_structures_from_geometry(
        rct_geos, prd_geos, zmat=zmat, stereo=stereo, ts_stereo=ts_stereo,
        ts_enant=ts_enant)


def with_structures_from_zmatrix(rct_zmas, prd_zmas, zmat=False, stereo=False,
                                 ts_stereo=None, ts_enant=None):
    """ Get reaction objects with geometry/z-matrix structures from SMILES

    :param rct_zmas: The reactant z-matrices
    :type rct_zmas: list of automol z-matrices
    :param prd_zmas: The product z-matrices
    :type prd_zmas: list of automol z-matrices
    :param zmat: Return z-matrix structures instead of cartesian geometries?
    :type zmat: bool
    :param stereo: Include stereoassignments? If `True`, then `ts_stereo` will
        default to `True` and `ts_enant` will default to `False`.
    :type stereo: bool
    :param ts_stereo: Include fleeting TS stereoisomers as distinct TSs?
    :type ts_stereo: bool or NoneType
    :param ts_enant: Include fleeting TS enantiomers as distinct TSs?
    :type ts_enant: bool or NoneType
    :returns: A series of tuples containing, in order, the reaction object, the
        TS structure, the (sorted) reactant structures, and the (sorted)
        product structures. The latter three are returned as either automol
        geometry objects or, if `zmatrix` was set to `True, automol z-matrix
        objects.
    """

    rct_geos = list(map(automol.zmat.geometry, rct_zmas))
    prd_geos = list(map(automol.zmat.geometry, prd_zmas))

    return with_structures_from_geometry(
        rct_geos, prd_geos, zmat=zmat, stereo=stereo, ts_stereo=ts_stereo,
        ts_enant=ts_enant)


def with_structures_from_geometry(rct_geos, prd_geos, zmat=False, stereo=False,
                                  ts_stereo=None, ts_enant=None):
    """ Get reaction objects with geometry/z-matrix structures from geometries

    :param rct_geos: The reactant geometries
    :type rct_geos: list of automol geometries
    :param prd_geos: The product geometries
    :type prd_geos: list of automol geometries
    :param zmat: Return z-matrix structures instead of cartesian geometries?
    :type zmat: bool
    :param stereo: Include stereoassignments? If `True`, then `ts_stereo` will
        default to `True` and `ts_enant` will default to `False`.
    :type stereo: bool
    :param ts_stereo: Include fleeting TS stereoisomers as distinct TSs?
    :type ts_stereo: bool or NoneType
    :param ts_enant: Include fleeting TS enantiomers as distinct TSs?
    :type ts_enant: bool or NoneType
    :returns: A series of tuples containing, in order, the reaction object, the
        TS structure, the (sorted) reactant structures, and the (sorted)
        product structures. The latter three are returned as either automol
        geometry objects or, if `zmatrix` was set to `True, automol z-matrix
        objects.
    """

    # Identify the reaction based on the reactants and products
    rct_gras = list(map(automol.geom.graph, rct_geos))
    prd_gras = list(map(automol.geom.graph, prd_geos))

    ts_enant = False if ts_enant is None else ts_enant
    if stereo:
        ts_stereo = True if ts_stereo is None else ts_stereo
    else:
        ts_stereo = False if ts_stereo is None else ts_stereo
        rct_gras = list(map(automol.graph.without_stereo_parities, rct_gras))
        prd_gras = list(map(automol.graph.without_stereo_parities, prd_gras))

    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = automol.reac.find(
        rct_gras, prd_gras, ts_stereo=ts_stereo, ts_enant=ts_enant)

    # Obtain the reaction objects and structures to return
    ret = tuple()
    for rxn in rxns:
        std_rxn, std_rgeos, std_pgeos = (
            automol.reac.standard_keys_with_sorted_geometries(
                rxn, rct_geos, prd_geos))

        # Add rxn object set to master list
        if std_rxn is not None:
            # Form the transition state geom using the rxn object
            ts_geo = automol.reac.ts_geometry(std_rxn, std_rgeos, log=False)

            # Determine which geometries to store
            if not zmat:
                ret += ((std_rxn, ts_geo, std_rgeos, std_pgeos),)
            else:
                ts_zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(
                    std_rxn, ts_geo)
                std_zrxn = automol.reac.relabel_for_zmatrix(
                    std_rxn, zma_keys, dummy_key_dct)
                rct_zmas = tuple(map(automol.geom.zmatrix, std_rgeos))
                prd_zmas = tuple(map(automol.geom.zmatrix, std_pgeos))

                ret += ((std_zrxn, ts_zma, rct_zmas, prd_zmas),)

    # Set to None if no objects found
    if not ret:
        ret = None

    return ret
