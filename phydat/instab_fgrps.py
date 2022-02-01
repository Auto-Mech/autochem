""" Dictionary that details what functional groups will yield
    electronically unstable species when they are attached to
    specified radical atoms. For a given X.-FGRP dissociation,
    we also denote one of the expected products.

DCT = {
    atm_symb: {
        instab_grp1: instab_prd1_from_grp1,
        instab_grp2: instab_prd1_from_grp2
    }
}
"""

DCT = {
    'C': {
        'InChI=1S/HO2/c1-2/h1H': 'InChI=1S/HO/h1H',
        'InChI=1S/NO3/c2-1(3)4': 'InChI=1S/NO2/c2-1-3'
    }
}
