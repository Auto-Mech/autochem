""" Storing a dictionary of active space sizes for various species
    DCT[(inchi, multiplicity)]
        = (num_active_orbs, num_active_elecs, num_states)
"""


DCT = {
    # Atoms
    ('InChI=1S/O', 3): (3, 4, 3),
    ('InChI=1S/O', 1): (3, 4, 5),
    ('InChI=1S/N', 3): (3, 3, 1),
    # Diatomics
    ('InChI=1S/O2/c1-2', 3): (4, 6, 1),
    ('InChI=1S/HO/h1H', 2): (2, 3, 2),
    # Molecules
    ('InChI=1S/CH2/h1H2', 3): (2, 2, 1),
    ('InChI=1S/CH2/h1H2', 1): (2, 2, 1)
}
