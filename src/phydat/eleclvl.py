"""
Library of electronic levels for several common, small species.
Dictionary formatted as:
    DCT[(InChI_Key, Multiplicity)] = [[energy(cm-1), degeneracy), ...)
"""

DCT = {
    ('InChI=1S/B', 2): ((0., 2), (16., 4)),
    ('InChI=1S/C', 3): ((0., 1), (16.4, 3), (43.5, 5)),
    ('InChI=1S/N', 2): ((0., 6), (8., 4)),
    ('InChI=1S/O', 3): ((0., 5), (158.5, 3), (226.5, 1)),
    ('InChI=1S/F', 2): ((0., 4), (404.1, 2)),
    ('InChI=1S/Cl', 2): ((0., 4), (883.4, 2)),
    ('InChI=1S/Br', 2): ((0., 4), (685.2, 2)),
    ('InChI=1S/HO/h1H', 2): ((0., 2), (138.9, 2)),
    ('InChI=1S/NO/c1-2', 2): ((0., 2), (123.1, 2)),
    ('InChI=1S/CH3O/c1-2/h1H3', 2): ((0., 3), (69.1, 1))
}
