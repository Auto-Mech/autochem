""" script to give a basic demo of distance geometry functionality
"""
import numpy
import automol


numpy.random.seed(2)

# ICH = ('InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)'
#        '4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/'
#        't10-,11+,13-,16-,17-/m0/s1')
# GEO = automol.inchi.geometry(ICH)
# GRA = automol.geom.graph(GEO)
GRA = automol.graph.from_string("""
    atoms:
      1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      3: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      4: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      5: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      6: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      7: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      8: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      9: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      10: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: true}
      11: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: true}
      12: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      13: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: true}
      14: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      15: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
      16: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: false}
      17: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: false}
      18: {symbol: N, implicit_hydrogen_valence: 0, stereo_parity: null}
      19: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
      20: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
      21: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
      22: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      23: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      24: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      25: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      26: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      27: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      28: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      29: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      30: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      31: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      32: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      33: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      34: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      35: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      36: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      37: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      38: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      39: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
      40: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
    bonds:
      1-18: {order: 1, stereo_parity: null}
      1-22: {order: 1, stereo_parity: null}
      1-23: {order: 1, stereo_parity: null}
      1-24: {order: 1, stereo_parity: null}
      2-4: {order: 1, stereo_parity: null}
      2-9: {order: 1, stereo_parity: null}
      2-25: {order: 1, stereo_parity: null}
      3-5: {order: 1, stereo_parity: null}
      3-10: {order: 1, stereo_parity: null}
      3-26: {order: 1, stereo_parity: null}
      4-12: {order: 1, stereo_parity: null}
      4-27: {order: 1, stereo_parity: null}
      5-13: {order: 1, stereo_parity: null}
      5-28: {order: 1, stereo_parity: null}
      6-7: {order: 1, stereo_parity: null}
      6-17: {order: 1, stereo_parity: null}
      6-29: {order: 1, stereo_parity: null}
      6-30: {order: 1, stereo_parity: null}
      7-18: {order: 1, stereo_parity: null}
      7-31: {order: 1, stereo_parity: null}
      7-32: {order: 1, stereo_parity: null}
      8-9: {order: 1, stereo_parity: null}
      8-11: {order: 1, stereo_parity: null}
      8-33: {order: 1, stereo_parity: null}
      8-34: {order: 1, stereo_parity: null}
      9-14: {order: 1, stereo_parity: null}
      10-11: {order: 1, stereo_parity: null}
      10-17: {order: 1, stereo_parity: null}
      10-35: {order: 1, stereo_parity: null}
      11-18: {order: 1, stereo_parity: null}
      11-36: {order: 1, stereo_parity: null}
      12-15: {order: 1, stereo_parity: null}
      12-19: {order: 1, stereo_parity: null}
      13-16: {order: 1, stereo_parity: null}
      13-20: {order: 1, stereo_parity: null}
      13-37: {order: 1, stereo_parity: null}
      14-15: {order: 1, stereo_parity: null}
      14-17: {order: 1, stereo_parity: null}
      15-21: {order: 1, stereo_parity: null}
      16-17: {order: 1, stereo_parity: null}
      16-21: {order: 1, stereo_parity: null}
      16-38: {order: 1, stereo_parity: null}
      19-39: {order: 1, stereo_parity: null}
      20-40: {order: 1, stereo_parity: null}
""")

KEYS = sorted(automol.graph.atom_keys(GRA))
LMAT, UMAT = automol.graph.embed.distance_bounds_matrices(GRA, KEYS)
XMAT = automol.graph.embed.sample_raw_distance_coordinates(GRA, KEYS,
                                                           dim4=True)

P_DCT = automol.graph.embed.planarity_constraint_bounds(GRA, KEYS)
XMAT, CONV = automol.embed.cleaned_up_coordinates(XMAT, LMAT, UMAT, P_DCT)
print(CONV)

SYMS = list(map(automol.graph.atom_symbols(GRA).__getitem__, KEYS))
GEO = automol.embed.geometry_from_coordinates(XMAT, SYMS)
print(automol.geom.string(GEO))

GRA = automol.graph.without_stereo_parities(GRA)
GRA2 = automol.geom.connectivity_graph(GEO)
print(GRA == GRA2)
