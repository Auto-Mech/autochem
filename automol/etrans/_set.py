""" Determine effective model
"""

import automol.inchi
import automol.geom
import automol.graph

# Break checks
BAD_ICHS = (
    'InChI=1S/H2/h1H'
)


def effective_model(well_ich, well_geo, bath_info):
    """ For the collision between a given well and bath species, determine
        which effective series would be the most suitable model for
        estimating the energy transfer parameters
    """

    # Initialize the models
    bath_model = bath_info[0]
    tgt_model = None
    # check if baths in series, if not set to Argon

    # Build the graph
    well_gra = automol.geom.graph(well_geo)

    # Identify the the target model
    if well_ich not in BAD_ICHS:
        # Set model based on broad values
        if automol.graph.radical_species(well_gra):
            tgt_model = '1-alkyl'
        elif automol.graph.hydrocarbon_species(well_gra):
            tgt_model = 'n-alkane'
        else:
            # Set priority based on bond-dissociation energies
            fgrp_dct = automol.graph.functional_group_dct(well_gra)
            if fgrp_dct[automol.graph.FunctionalGroup.HYDROPEROXY]:
                tgt_model = 'n-hydroperoxide'
            elif fgrp_dct[automol.graph.FunctionalGroup.EPOXIDE]:
                tgt_model = 'epoxide'
            elif fgrp_dct[automol.graph.FunctionalGroup.ETHER]:
                tgt_model = 'ether'
            elif fgrp_dct[automol.graph.FunctionalGroup.ALCOHOL]:
                tgt_model = 'n-alcohol'

        # For now, set model to alkanes if nothing found and set up return obj
        if tgt_model is None:
            tgt_model = 'n-alkane'

        ret = (bath_model, tgt_model)
    else:
        ret = None

    return ret
