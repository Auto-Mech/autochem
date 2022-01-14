""" testing fake_stereo_geometry function
"""

import automol

# ICH = automol.smiles.inchi('FC=C(F)N')
# ICH = automol.smiles.inchi('FC=CC=CCCCC(O)(C)')
# ICH = ('InChI=1S/C7H14O3/c1-3-4-6-7(10-8)5(2)9-6/h5-8H,3-4H2,1-2H3/'
#        't5-,6-,7-/m0/s1')
# ICH = ('InChI=1S/C7H14O3/c1-3-6(10-8)7-4-5(2)9-7/h5-8H,3-4H2,1-2H3/'
#        't5-,6+,7+/m1/s1')
# ICH = ('InChI=1S/C8H13/c1-5-7(3)8(4)6-2/h5-7H,1-2H2,3-4H3/t7-/m0/s1')
# ICH = ('InChI=1S/C8H13/c1-3-5-7-8-6-4-2/h3-7H,8H2,1-2H3/b6-4-')
ICH = ('InChI=1S/C7H14O3/c1-5-3-7(9-5)4-6(2)10-8/h5-8H,3-4H2,1-2H3/'
       't5-,6+,7-/m0/s1')
GEO_IN = automol.inchi.geometry(ICH)
ICH_REF = automol.geom.inchi(GEO_IN)
print(ICH_REF)
GRA = automol.geom.graph(GEO_IN)
ICH = automol.graph.inchi(GRA)
print(ICH)
assert ICH_REF == ICH
