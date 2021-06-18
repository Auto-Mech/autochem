""" Determine effective model
"""

import automol.inchi
import automol.geom
import automol.graph
from automol.etrans._par import D0_GRP_LST


# Break checks
BAD_ICHS = (
    'InChI=1S/H2/h1H'
)


def effective_model(tgt_ich, bath_ich):
    """ For the collision between a given tgt and bath species, determine
        which effective series would be the most suitable model for
        estimating the energy transfer parameters
    """

    # Initialize the model
    tgt_model = None

    # Build the graph
    tgt_gra = automol.geom.graph(automol.inchi.geometry(tgt_ich))

    # Identify the the target model
    if tgt_ich not in BAD_ICHS:
        # Set model based on broad values
        if automol.graph.radical_species(tgt_gra):
            tgt_model = '1-alkyl'
        elif automol.graph.hydrocarbon_species(tgt_gra):
            tgt_model = 'n-alkane'
        else:
            # Set priority based on bond-dissociation energies
            fgrp_dct = automol.graph.functional_group_dct(tgt_gra)
            for (fgrp, model) in D0_GRP_LST:
                if fgrp_dct[fgrp]:
                    tgt_model = model
                    break

        # For now, set model to alkanes if nothing found and set up return obj
        if tgt_model is None:
            tgt_model = 'n-alkane'

        ret = (bath_ich, tgt_model)
    else:
        ret = None

    return ret
