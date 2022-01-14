"""
  Test geom comparison functions
"""


import automol.geom


GEO1 = (
    ('C', (-1.4283035320563338, 0.013425343735546437, -0.030302158896694683)),
    ('C', (1.4283027358735494, -0.013425597530894248, 0.0303022919384165)),
    ('H', (-2.1972722614281355, -0.19229727219177065, 1.8778380427620682)),
    ('H', (-2.121310184939721, 1.792702413487708, -0.8231106338374065)),
    ('H', (-2.1448124562913287, -1.5396513482615042, -1.191852168914227)),
    ('H', (2.1448121742707795, 1.539654946791746, 1.1918517388178247)),
    ('H', (2.1972712765396953, 0.1922944277301287, -1.8778395029874426)),
    ('H', (2.121312248031497, -1.7927029137609576, 0.8231123911174519)))
# geom that is similar to GEO1
GEO1_DISP = (
    ('C', (8.571696467943667, 10.013425343735546, 9.969697841103306)),
    ('C', (11.42830273587355, 9.986574402469106, 10.030302291938417)),
    ('H', (7.802727738571864, 9.807702727808229, 11.877838042762068)),
    ('H', (7.878689815060278, 11.792702413487708, 9.176889366162593)),
    ('H', (7.855187543708672, 8.460348651738496, 8.808147831085773)),
    ('H', (12.144812174270779, 11.539654946791746, 11.191851738817824)),
    ('H', (12.197271276539695, 10.192294427730129, 8.122160497012558)),
    ('H', (12.121312248031497, 8.207297086239043, 10.823112391117451)))
# geoms that are NOT similar to GEO1
GEO1_DIFF_ATOM_ORDER = (
    ('C', (-1.4283035320563338, 0.013425343735546437, -0.030302158896694683)),
    ('H', (-2.1972722614281355, -0.19229727219177065, 1.8778380427620682)),
    ('H', (-2.121310184939721, 1.792702413487708, -0.8231106338374065)),
    ('H', (-2.1448124562913287, -1.5396513482615042, -1.191852168914227)),
    ('C', (1.4283027358735494, -0.013425597530894248, 0.0303022919384165)),
    ('H', (2.1448121742707795, 1.539654946791746, 1.1918517388178247)),
    ('H', (2.1972712765396953, 0.1922944277301287, -1.8778395029874426)),
    ('H', (2.121312248031497, -1.7927029137609576, 0.8231123911174519)))
GEO1_CC_STRETCH = (
    ('C', (-2.5382105317, 0.0134246134, -0.0303017562)),
    ('C', (1.4283021453, -0.0134265031, 0.0303017562)),
    ('H', (-3.3071778309, -0.1922966269, 1.8778377224)),
    ('H', (-3.2312165154, 1.7927017871, -0.8231098386)),
    ('H', (-3.2547190375, -1.5396504697, -1.1918520707)),
    ('H', (2.1448125409, 1.5396542492, 1.1918520707)),
    ('H', (2.1972713343, 0.1922947371, -1.8778396121)),
    ('H', (2.1213119085, -1.7927036769, 0.8231117283)))

GEO_LST1 = (GEO1, GEO1_DISP)
GEO_LST2 = (GEO1_DIFF_ATOM_ORDER, GEO1_DISP)
GEO_LST3 = (GEO1_DIFF_ATOM_ORDER,)
GEO_LST4 = (GEO1_CC_STRETCH,)
GEO_LST5 = (GEO1_DIFF_ATOM_ORDER, GEO1_CC_STRETCH)

CHECK_DCT = {
    'dist': None, 'stereo': None, 'coulomb': None}
CHECK_DCT2 = {
    'dist': 3.5e-1, 'stereo': None, 'coulomb': 1.5e-2}


# def test__newzma():
#     """ test
#     """
#
#     build_remdummy_shift_lst(zma)
#     shift_vals_from_dummy(vals, zma)
#     is_atom_closest_to_bond_atom(zma, idx_rad, bond_dist)
#     calc_rxn_angle(ts_zma, frm_bnd_keys, brk_bnd_keys, rxn_class)
#
#
def test__comp():
    """ test automol.comp.is_unique
    """

    # Test CHECK DCT with values
    unique1, idx1 = automol.geom.is_unique(GEO1, GEO_LST1, CHECK_DCT)
    unique2, idx2 = automol.geom.is_unique(GEO1, GEO_LST2, CHECK_DCT)
    unique3, idx3 = automol.geom.is_unique(GEO1, GEO_LST3, CHECK_DCT)
    unique4, idx4 = automol.geom.is_unique(GEO1, GEO_LST4, CHECK_DCT)
    unique5, idx5 = automol.geom.is_unique(GEO1, GEO_LST5, CHECK_DCT)
    assert not unique1 and idx1 == 0
    assert not unique2 and idx2 == 1
    assert unique3 and idx3 is None
    assert unique4 and idx4 is None
    assert unique5 and idx5 is None

    # Test CHECK DCT with values
    unique6, idx6 = automol.geom.is_unique(GEO1, GEO_LST1, CHECK_DCT2)
    assert not unique6 and idx6 == 0


if __name__ == '__main__':
    test__comp()
