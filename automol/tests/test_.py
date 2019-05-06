""" test the automechanic.mol module
"""
import os
import numpy
import automol as mol

PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')

HEPTANE_ICHS = numpy.loadtxt(os.path.join(DATA_PATH, 'heptane_inchis.txt'),
                             dtype=str)

AR_ICH = 'InChI=1S/Ar'

C2H2F2_SMI = 'F/C=C/F'
C2H2F2_SMI_NO_STEREO = 'FC=CF'
C2H2F2_ICH = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+'
C2H2F2_ICH_NO_STEREO = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H'
C2H2F2_ICH_STEREO_UNKNOWN = 'InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1?'
C2H2F2_GEO = (('F', (2.994881276150, -1.414434615111, -0.807144415388)),
              ('C', (1.170155936996, 0.359360756989, -0.513323178859)),
              ('C', (-1.201356763194, -0.347546894407, -0.3408392500119)),
              ('F', (-3.027970874978, 1.39211904938, -0.0492290974807)),
              ('H', (1.731596406235, 2.324260256203, -0.4292070203467)),
              ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))

C8H13O_ICH = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-/m0/s1')
C8H13O_ICH_NO_ENANTIOMER = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-')
C8H13O_ICH_PARTIAL_STEREO = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/b5-3-/t8-/m0/s1')
C8H13O_ICH_NO_STEREO = 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3'
C8H13O_ICHS = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-/m0/s1',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4+/t8-/m0/s1',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3+,6-4-/t8-/m0/s1',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3+,6-4+/t8-/m0/s1'
)

C2H4CLF_SMI_NO_STEREO = 'C(Cl)(F)C'
C2H4CLF_ICH = 'InChI=1S/C2H4ClF/c1-2(3)4/h2H,1H3/t2-/m0/s1'
C2H4CLF_ICH_NO_STEREO = 'InChI=1S/C2H4ClF/c1-2(3)4/h2H,1H3'
C2H4CLF_ICH_INCOMPLETE_STEREO = 'InChI=1S/C2H4ClF/c1-2(3)4/h2H,1H3/t2-'
C2H4CLF_ICH_STEREO_UNKNOWN = 'InChI=1/C2H4ClF/c1-2(3)4/h2H,1H3/t2?'

# rdkit geometry failures:
RDKIT_FAIL_ICHS = (
    'InChI=1S/C6H10/c1-3-5-6-4-2/h3-6H,1-2H3',
    'InChI=1S/C5H10/c1-3-5-4-2/h3,5H,4H2,1-2H3',
    'InChI=1S/C4H7O2/c1-4(2)3-6-5/h1,3H2,2H3',
    'InChI=1S/C8H13O/c1-2-3-4-5-6-7-8-9/h2-3,6-7H,4-5,8H2,1H3')
# pybel geometry failures:
PYBEL_FAIL_ICHS = (
    'InChI=1S/C3H7O/c1-3(2)4/h3H,1-2H3',
    'InChI=1S/C3H7O2/c1-3(2)5-4/h3-4H,1H2,2H3',
    'InChI=1S/C3H7O4/c4-6-2-1-3-7-5/h4H,1-3H2',
    'InChI=1S/C3H6O3/c4-2-1-3-6-5/h2,5H,1,3H2',
    'InChI=1S/C3H7O4/c4-6-2-1-3-7-5/h1,4-5H,2-3H2',
    'InChI=1S/C3H6O3/c4-6-3-1-5-2-3/h3-4H,1-2H2',
    'InChI=1S/C3H6O/c1-2-4-3-1/h1-3H2',
    'InChI=1S/C3H2/c1-2-3-1/h1-2H',
    'InChI=1S/C4H9O/c1-3-4(2)5/h4H,3H2,1-2H3',
    'InChI=1S/C4H9O2/c1-3-4(2)6-5/h4H,3H2,1-2H3',
    'InChI=1S/C4H10O2/c1-3-4(2)6-5/h4-5H,3H2,1-2H3',
    'InChI=1S/C4H8O/c1-4-2-3-5-4/h4H,2-3H2,1H3',
    'InChI=1S/C4H8O/c1-2-4-5-3-1/h1-4H2',
    'InChI=1S/C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3',
    'InChI=1S/C4H9O4/c1-4(8-6)2-3-7-5/h4-5H,2-3H2,1H3',
    'InChI=1S/C4H9O4/c5-7-3-1-2-4-8-6/h5H,1-4H2',
    'InChI=1S/C4H9O4/c1-3(7-5)4(2)8-6/h3-6H,1H2,2H3',
    'InChI=1S/C4H8O3/c1-3-4(7-5)2-6-3/h3-5H,2H2,1H3',
    'InChI=1S/C4H8O3/c5-7-2-1-4-3-6-4/h4-5H,1-3H2',
    'InChI=1S/C4H8O3/c5-7-4-1-2-6-3-4/h4-5H,1-3H2',
    'InChI=1S/C4H8O3/c1-3(7-5)4-2-6-4/h3-5H,2H2,1H3',
    'InChI=1S/C4H8O3/c1-3-4(7-3)2-6-5/h3-5H,2H2,1H3',
    'InChI=1S/C4H8O3/c1-4(7-6)2-3-5/h3-4,6H,2H2,1H3',
    'InChI=1S/C4H8O3/c5-3-1-2-4-7-6/h3,6H,1-2,4H2',
    'InChI=1S/C4H9O/c1-4(2)3-5/h4H,3H2,1-2H3',
    'InChI=1S/C4H9O2/c1-4(2,3)6-5/h1-3H3',
    'InChI=1S/C4H10O2/c1-4(2,3)6-5/h5H,1-3H3',
    'InChI=1S/C4H8O4/c1-4(8-6)2-7-3(4)5/h3,5-6H,2H2,1H3',
    'InChI=1S/C4H9O3/c1-2-4(5)3-7-6/h4,6H,2-3H2,1H3',
    'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3',
    'InChI=1S/C4H7O3/c5-3-1-2-4-7-6/h1-2,6H,3-4H2',
    'InChI=1S/C4H7O2/c1-2-3-4-6-5/h3-4H,2H2,1H3',
    'InChI=1S/C4H9O3/c1-2-4(5)3-7-6/h4-5H,2-3H2,1H3',
    'InChI=1S/C4H6O/c1-2-4-3-5-4/h2,4H,1,3H2',
    'InChI=1S/C5H12O2/c1-3-4-5(2)7-6/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C5H12O2/c1-3-5(4-2)7-6/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C5H11O2/c1-3-4-5(2)7-6/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O2/c1-3-5(4-2)7-6/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O/c1-3-4-5(2)6/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O/c1-3-5(6)4-2/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O2/c1-3-5(4-2)7-6/h5-6H,1,3-4H2,2H3',
    'InChI=1S/C5H10O/c1-2-5-3-4-6-5/h5H,2-4H2,1H3',
    'InChI=1S/C5H10O/c1-5-3-2-4-6-5/h5H,2-4H2,1H3',
    'InChI=1S/C5H10O/c1-2-4-6-5-3-1/h1-5H2',
    'InChI=1S/C5H10O/c1-3-5-4(2)6-5/h4-5H,3H2,1-2H3',
    'InChI=1S/C5H10O/c1-4-3-5(2)6-4/h4-5H,3H2,1-2H3',
    'InChI=1S/C5H11O4/c1-2-5(9-7)3-4-8-6/h5-6H,2-4H2,1H3',
    'InChI=1S/C5H11O4/c1-5(9-7)3-2-4-8-6/h5-6H,2-4H2,1H3',
    'InChI=1S/C5H11O4/c6-8-4-2-1-3-5-9-7/h6H,1-5H2',
    'InChI=1S/C5H11O4/c1-4(8-6)3-5(2)9-7/h4-6H,3H2,1-2H3',
    'InChI=1S/C5H11O4/c1-2-5(9-7)3-4-8-6/h5,7H,2-4H2,1H3',
    'InChI=1S/C5H11O4/c1-3-5(9-7)4(2)8-6/h4-5,7H,3H2,1-2H3',
    'InChI=1S/C5H11O4/c1-5(9-7)3-2-4-8-6/h5-7H,1-4H2',
    'InChI=1S/C5H11O4/c1-3-5(9-7)4(2)8-6/h4-7H,1,3H2,2H3',
    'InChI=1S/C5H11O4/c1-2-5(9-7)3-4-8-6/h5-7H,1-4H2',
    'InChI=1S/C5H11O4/c1-2-3-5(9-7)4-8-6/h5-7H,1-4H2',
    'InChI=1S/C5H11O4/c6-8-4-2-1-3-5-9-7/h2,6-7H,1,3-5H2',
    'InChI=1S/C5H11O4/c1-2-3-5(9-7)4-8-6/h2,5-7H,3-4H2,1H3',
    'InChI=1S/C5H11O4/c1-5(9-7)3-2-4-8-6/h3,5-7H,2,4H2,1H3',
    'InChI=1S/C5H10O2/c1-3-4-5(2)7-6/h3,5-6H,1,4H2,2H3',
    'InChI=1S/C5H10O2/c1-3-4-5(2)7-6/h3-6H,1-2H3',
    'InChI=1S/C5H10O3/c1-2-4-5(8-6)3-7-4/h4-6H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c1-4(8-6)5-2-3-7-5/h4-6H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c1-4-2-5(8-6)3-7-4/h4-6H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c6-8-4-5-2-1-3-7-5/h5-6H,1-4H2',
    'InChI=1S/C5H10O3/c6-8-5-2-1-3-7-4-5/h5-6H,1-4H2',
    'InChI=1S/C5H10O3/c1-3-5(7-3)4(2)8-6/h3-6H,1-2H3',
    'InChI=1S/C5H10O3/c1-4-5(8-4)2-3-7-6/h4-6H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c1-2-5(8-7)3-4-6/h4-5,7H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c1-5(8-7)3-2-4-6/h4-5,7H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c6-4-2-1-3-5-8-7/h4,7H,1-3,5H2',
    'InChI=1S/C5H10O3/c1-4(6)3-5(2)8-7/h5,7H,3H2,1-2H3',
    'InChI=1S/C5H10O3/c1-2-5(6)3-4-8-7/h7H,2-4H2,1H3',
    'InChI=1S/C5H10O3/c1-3-5(6)4(2)8-7/h4,7H,3H2,1-2H3',
    'InChI=1S/C5H9O2/c1-2-5(7)3-4-6/h4-5H,2-3H2,1H3',
    'InChI=1S/C5H9O2/c1-5(7)3-2-4-6/h4-5H,2-3H2,1H3',
    'InChI=1S/C5H9O2/c6-4-2-1-3-5-7/h4H,1-3,5H2',
    'InChI=1S/C5H9O2/c1-5(7)3-2-4-6/h2-4H2,1H3',
    'InChI=1S/C5H9O2/c1-2-5(7)3-4-6/h2-4H2,1H3',
    'InChI=1S/C5H9O2/c1-3-5(7)4(2)6/h4H,3H2,1-2H3',
    'InChI=1S/C5H8O2/c1-4(6)3-5(2)7/h3H2,1-2H3',
    'InChI=1S/C5H12O3/c1-4(6)3-5(2)8-7/h4-7H,3H2,1-2H3',
    'InChI=1S/C5H12O3/c1-2-5(8-7)3-4-6/h5-7H,2-4H2,1H3',
    'InChI=1S/C5H10O3/c1-4-3-5(2,6)8-7-4/h4,6H,3H2,1-2H3',
    'InChI=1S/C5H10O3/c1-2-5(6)3-4-7-8-5/h6H,2-4H2,1H3',
    'InChI=1S/C5H12O2/c1-4-5(2,3)7-6/h6H,4H2,1-3H3',
    'InChI=1S/C5H12O2/c1-4(2)5(3)7-6/h4-6H,1-3H3',
    'InChI=1S/C5H12O2/c1-5(2)3-4-7-6/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C5H11O2/c1-4-5(2,3)7-6/h4H2,1-3H3',
    'InChI=1S/C5H11O2/c1-4(2)5(3)7-6/h4-5H,1-3H3',
    'InChI=1S/C5H11O2/c1-5(2)3-4-7-6/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O/c1-4-5(2,3)6/h4H2,1-3H3',
    'InChI=1S/C5H11O/c1-4(2)5(3)6/h4-5H,1-3H3',
    'InChI=1S/C5H11O/c1-5(2)3-4-6/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O2/c1-4(2)5(3)7-6/h4-6H,1H2,2-3H3',
    'InChI=1S/C5H10O/c1-5-2-3-6-4-5/h5H,2-4H2,1H3',
    'InChI=1S/C5H10O/c1-5(2)3-4-6-5/h3-4H2,1-2H3',
    'InChI=1S/C5H11O4/c1-4(8-6)5(2,3)9-7/h4,7H,1-3H3',
    'InChI=1S/C5H11O4/c1-5(2,9-7)3-4-8-6/h7H,3-4H2,1-2H3',
    'InChI=1S/C5H11O4/c1-5(4-9-7)2-3-8-6/h5-7H,1-4H2',
    'InChI=1S/C5H10O3/c6-8-2-1-5-3-7-4-5/h5-6H,1-4H2',
    'InChI=1S/C5H10O3/c1-4(8-6)5(2)3-7-5/h4,6H,3H2,1-2H3',
    'InChI=1S/C5H10O3/c1-5(4-7-5)2-3-8-6/h6H,2-4H2,1H3',
    'InChI=1S/C5H10O3/c6-8-4-5-1-2-7-3-5/h5-6H,1-4H2',
    'InChI=1S/C5H10O3/c1-5(2)4(8-5)3-7-6/h4,6H,3H2,1-2H3',
    'InChI=1S/C5H10O3/c1-5(4-8-6)2-3-7-5/h6H,2-4H2,1H3',
    'InChI=1S/C5H10O3/c1-5(2)4(8-6)3-7-5/h4,6H,3H2,1-2H3',
    'InChI=1S/C5H11O3/c1-4(8-7)5(2,3)6/h4,6H,1-3H3',
    'InChI=1S/C5H9O/c1-4-5(2,3)6/h4H,1H2,2-3H3',
    'InChI=1S/C5H9O3/c1-2-5(8-7)3-4-6/h4-5H,2-3H2,1H3',
    'InChI=1S/C5H9O3/c1-5(8-7)3-2-4-6/h4-5H,2-3H2,1H3',
    'InChI=1S/C5H9O3/c6-4-2-1-3-5-8-7/h4H,1-3,5H2',
    'InChI=1S/C5H10O/c1-3-5(6)4-2/h3-4H2,1-2H3',
    'InChI=1S/C5H9O/c1-3-5(6)4-2/h1,3-4H2,2H3',
    'InChI=1S/C5H9O3/c1-4(6)3-5(2)8-7/h5H,3H2,1-2H3',
    'InChI=1S/C5H9O3/c1-2-5(6)3-4-8-7/h2-4H2,1H3',
    'InChI=1S/C5H9O3/c1-3-5(6)4(2)8-7/h4H,3H2,1-2H3',
    'InChI=1S/C5H8O/c1-3-5(6)4-2/h3H,1,4H2,2H3',
    'InChI=1S/C5H11O3/c1-2-3-5(4-6)8-7/h5-6H,2-4H2,1H3',
    'InChI=1S/C5H11O3/c1-2-5(8-7)3-4-6/h5-6H,2-4H2,1H3',
    'InChI=1S/C5H11O3/c1-5(8-7)3-2-4-6/h5-6H,2-4H2,1H3',
    'InChI=1S/C5H11O3/c6-4-2-1-3-5-8-7/h6H,1-5H2',
    'InChI=1S/C6H14O2/c1-3-4-5-6(2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C6H14O2/c1-3-5-6(4-2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C6H13O2/c1-3-4-5-6(2)8-7/h6H,3-5H2,1-2H3',
    'InChI=1S/C6H13O2/c1-3-5-6(4-2)8-7/h6H,3-5H2,1-2H3',
    'InChI=1S/C6H13O/c1-3-4-5-6(2)7/h6H,3-5H2,1-2H3',
    'InChI=1S/C6H13O/c1-3-5-6(7)4-2/h6H,3-5H2,1-2H3',
    'InChI=1S/C6H13O2/c1-3-4-5-6(2)8-7/h6-7H,2-5H2,1H3',
    'InChI=1S/C6H12O/c1-2-3-6-4-5-7-6/h6H,2-5H2,1H3',
    'InChI=1S/C6H12O/c1-2-6-4-3-5-7-6/h6H,2-5H2,1H3',
    'InChI=1S/C6H12O/c1-6-4-2-3-5-7-6/h6H,2-5H2,1H3',
    'InChI=1S/C6H12O/c1-3-4-6-5(2)7-6/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C6H12O/c1-3-6-4-5(2)7-6/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C6H12O/c1-5-3-4-6(2)7-5/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C6H12O/c1-3-5-6(4-2)7-5/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C6H13O4/c1-2-3-6(10-8)4-5-9-7/h6-7H,2-5H2,1H3',
    'InChI=1S/C6H13O4/c1-2-6(10-8)4-3-5-9-7/h6-7H,2-5H2,1H3',
    'InChI=1S/C6H13O4/c1-6(10-8)4-2-3-5-9-7/h6-7H,2-5H2,1H3',
    'InChI=1S/C6H13O4/c1-3-6(10-8)4-5(2)9-7/h5-7H,3-4H2,1-2H3',
    'InChI=1S/C6H13O4/c1-5(9-7)3-4-6(2)10-8/h5-7H,3-4H2,1-2H3',
    'InChI=1S/C6H13O4/c1-2-3-6(10-8)4-5-9-7/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O4/c1-3-4-6(10-8)5(2)9-7/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O4/c1-2-6(10-8)4-3-5-9-7/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H12O2/c1-3-4-5-6(2)8-7/h3,6-7H,1,4-5H2,2H3',
    'InChI=1S/C6H13O4/c1-5(9-7)3-4-6(2)10-8/h5-8H,1,3-4H2,2H3',
    'InChI=1S/C6H12O3/c1-2-3-6(9-8)4-5-7/h5-6,8H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-2-6(9-8)4-3-5-7/h5-6,8H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-6(9-8)4-2-3-5-7/h5-6,8H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-3-6(9-8)4-5(2)7/h6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H12O3/c1-5(7)3-4-6(2)9-8/h6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H12O3/c1-2-3-6(7)4-5-9-8/h8H,2-5H2,1H3',
    'InChI=1S/C6H12O3/c1-3-4-6(7)5(2)9-8/h5,8H,3-4H2,1-2H3',
    'InChI=1S/C6H12O3/c1-2-6(7)4-3-5-9-8/h8H,2-5H2,1H3',
    'InChI=1S/C6H12O3/c1-5(9-7)2-3-6-4-8-6/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-2-3-5-6(9-7)4-8-5/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-2-5(9-7)6-3-4-8-6/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-2-5-3-6(9-7)4-8-5/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-5(9-7)6-3-2-4-8-6/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-5-2-3-6(9-7)4-8-5/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c7-9-5-6-3-1-2-4-8-6/h6-7H,1-5H2',
    'InChI=1S/C6H12O3/c1-2-3-5-6(9-5)4-8-7/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-3-5(9-7)6-4(2)8-6/h4-7H,3H2,1-2H3',
    'InChI=1S/C6H12O3/c1-4(9-7)3-6-5(2)8-6/h4-7H,3H2,1-2H3',
    'InChI=1S/C6H12O3/c1-5-6(9-5)3-2-4-8-7/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-4-3-6(8-4)5(2)9-7/h4-7H,3H2,1-2H3',
    'InChI=1S/C6H12O3/c1-2-5-6(9-5)3-4-8-7/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-3-5-6(8-5)4(2)9-7/h4-7H,3H2,1-2H3',
    'InChI=1S/C4H7O/c1-4-2-3-5-4/h4H,1-3H2',
    'InChI=1S/C4H7O/c1-3-4(2)5-3/h3-4H,1H2,2H3',
    'InChI=1S/C5H9O/c1-2-5-3-4-6-5/h5H,1-4H2',
    'InChI=1S/C5H9O/c1-5-3-2-4-6-5/h5H,1-4H2',
    'InChI=1S/C5H9O/c1-4-3-5(2)6-4/h4-5H,1,3H2,2H3',
    'InChI=1S/C6H11O2/c1-2-3-6(8)4-5-7/h5-6H,2-4H2,1H3',
    'InChI=1S/C6H11O2/c1-2-6(8)4-3-5-7/h5-6H,2-4H2,1H3',
    'InChI=1S/C6H11O2/c1-6(8)4-2-3-5-7/h5-6H,2-4H2,1H3',
    'InChI=1S/C6H11O2/c1-3-6(8)4-5(2)7/h6H,3-4H2,1-2H3',
    'InChI=1S/C6H11O2/c1-6(8)4-2-3-5-7/h2-5H2,1H3',
    'InChI=1S/C6H11O2/c1-2-3-6(8)4-5-7/h2-5H2,1H3',
    'InChI=1S/C6H11O2/c1-3-4-6(8)5(2)7/h5H,3-4H2,1-2H3',
    'InChI=1S/C6H11O2/c1-3-5(7)6(8)4-2/h5H,3-4H2,1-2H3',
    'InChI=1S/C6H11O2/c1-3-6(8)4-5(2)7/h5H,3-4H2,1-2H3',
    'InChI=1S/C6H11O2/c1-2-6(8)4-3-5-7/h2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-2-3-6(7)4-5-9-8/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-2-6(7)4-3-5-9-8/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-6(7)4-2-3-5-9-8/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-2-3-4-6(5-7)9-8/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-3-4-6(7)5(2)9-8/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-3-6(7)4-5(2)9-8/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-5(7)3-4-6(2)9-8/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-6(9-8)4-2-3-5-7/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-2-3-6(9-8)4-5-7/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-3-4-6(9-8)5(2)7/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-3-5(7)6(4-2)9-8/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-3-6(9-8)4-5(2)7/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-2-6(9-8)4-3-5-7/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H12O/c1-3-4-5-6(2)7/h3-5H2,1-2H3',
    'InChI=1S/C6H13O3/c1-2-3-4-6(7)5-9-8/h6-7H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-3-4-6(9-8)5(2)7/h5-7H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-3-4-6(7)5(2)9-8/h5-7H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-3-5(7)6(4-2)9-8/h5-7H,3-4H2,1-2H3',
    'InChI=1S/C6H11O/c1-3-4-5-6(2)7/h3,6H,1,4-5H2,2H3',
    'InChI=1S/C7H16O2/c1-3-4-5-6-7(2)9-8/h7-8H,3-6H2,1-2H3',
    'InChI=1S/C7H16O2/c1-3-5-6-7(4-2)9-8/h7-8H,3-6H2,1-2H3',
    'InChI=1S/C7H16O2/c1-3-5-7(9-8)6-4-2/h7-8H,3-6H2,1-2H3',
    'InChI=1S/C7H15O2/c1-3-4-5-6-7(2)9-8/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O2/c1-3-5-6-7(4-2)9-8/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O2/c1-3-5-7(9-8)6-4-2/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O/c1-3-4-5-6-7(2)8/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O/c1-3-5-6-7(8)4-2/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O/c1-3-5-7(8)6-4-2/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O2/c1-3-5-7(9-8)6-4-2/h7-8H,1,3-6H2,2H3',
    'InChI=1S/C7H14O/c1-2-3-4-7-5-6-8-7/h7H,2-6H2,1H3',
    'InChI=1S/C7H14O/c1-2-4-7-5-3-6-8-7/h7H,2-6H2,1H3',
    'InChI=1S/C7H14O/c1-2-7-5-3-4-6-8-7/h7H,2-6H2,1H3',
    'InChI=1S/C7H14O/c1-3-4-5-7-6(2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H14O/c1-3-4-7-5-6(2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H14O/c1-3-7-5-4-6(2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H14O/c1-6-4-3-5-7(2)8-6/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H14O/c1-3-5-7-6(4-2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H14O/c1-3-6-5-7(4-2)8-6/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H15O4/c1-2-3-4-7(11-9)5-6-10-8/h7-8H,2-6H2,1H3',
    'InChI=1S/C7H15O4/c1-2-4-7(11-9)5-3-6-10-8/h7-8H,2-6H2,1H3')


def test__smiles__inchi():
    """ smiles.inchi
    """
    assert mol.smiles.inchi(C2H2F2_SMI) == C2H2F2_ICH


def test__inchi__smiles():
    """ inchi.smiles
    """
    assert mol.smiles.inchi(mol.inchi.smiles(AR_ICH)) == AR_ICH
    assert (mol.smiles.inchi(mol.inchi.smiles(C8H13O_ICH_NO_STEREO))
            == C8H13O_ICH_NO_STEREO)
    assert (tuple(map(mol.smiles.inchi, map(mol.inchi.smiles, C8H13O_ICHS)))
            == C8H13O_ICHS)


def test__inchi__recalculate():
    """ inchi.recalculate
    """
    assert mol.inchi.recalculate(C2H2F2_ICH_NO_STEREO) == C2H2F2_ICH_NO_STEREO
    assert (mol.inchi.recalculate(C2H2F2_ICH_NO_STEREO, force_stereo=True)
            == C2H2F2_ICH_STEREO_UNKNOWN)


def test__inchi__is_closed():
    """ inchi.is_closed
    """
    assert mol.inchi.is_closed(C8H13O_ICH) is True
    assert mol.inchi.is_closed(C8H13O_ICH_PARTIAL_STEREO) is False
    assert mol.inchi.is_closed(C8H13O_ICH_NO_STEREO) is True
    assert mol.inchi.is_closed(C8H13O_ICH_NO_ENANTIOMER) is False


def test__inchi__prefix():
    """ inchi.prefix
    """
    assert mol.inchi.prefix(C2H2F2_ICH) == 'InChI=1S'
    assert mol.inchi.prefix(C2H2F2_ICH_NO_STEREO) == 'InChI=1S'
    assert mol.inchi.prefix(C2H2F2_ICH_STEREO_UNKNOWN) == 'InChI=1'


def test__inchi__version():
    """ inchi.version
    """
    assert mol.inchi.version(C2H2F2_ICH) == '1S'
    assert mol.inchi.version(C2H2F2_ICH_NO_STEREO) == '1S'
    assert mol.inchi.version(C2H2F2_ICH_STEREO_UNKNOWN) == '1'


def test__inchi__formula_layer():
    """ inchi.formula_layer
    """
    assert mol.inchi.formula_layer(C2H2F2_ICH) == 'C2H2F2'
    assert (mol.inchi.formula_layer('InChI=1S/2C2H5.Zn/c2*1-2;/h2*1H2,2H3;')
            == '2C2H5.Zn')


def test__inchi__key_layer():
    """ inchi.key_layer
    """
    assert mol.inchi.key_layer(C2H2F2_ICH, 'c') == 'c3-1-2-4'
    assert mol.inchi.key_layer(C2H2F2_ICH, 'h') == 'h1-2H'
    assert mol.inchi.key_layer(C2H2F2_ICH, 'b') == 'b2-1+'
    assert mol.inchi.key_layer(C2H2F2_ICH_STEREO_UNKNOWN, 'c') == 'c3-1-2-4'
    assert mol.inchi.key_layer(C2H2F2_ICH_STEREO_UNKNOWN, 'h') == 'h1-2H'
    assert mol.inchi.key_layer(C2H2F2_ICH_STEREO_UNKNOWN, 'b') == 'b2-1?'
    assert mol.inchi.key_layer(C2H2F2_ICH_NO_STEREO, 'c') == 'c3-1-2-4'
    assert mol.inchi.key_layer(C2H2F2_ICH_NO_STEREO, 'h') == 'h1-2H'
    assert mol.inchi.key_layer(C2H2F2_ICH_NO_STEREO, 'b') is None


def test__inchi__key_layer_content():
    """ inchi.key_layer_content
    """
    assert mol.inchi.key_layer_content(C2H2F2_ICH, 'c') == '3-1-2-4'
    assert mol.inchi.key_layer_content(C2H2F2_ICH, 'h') == '1-2H'
    assert mol.inchi.key_layer_content(C2H2F2_ICH, 'b') == '2-1+'
    assert (mol.inchi.key_layer_content(C2H2F2_ICH_STEREO_UNKNOWN, 'c')
            == '3-1-2-4')
    assert (mol.inchi.key_layer_content(C2H2F2_ICH_STEREO_UNKNOWN, 'h')
            == '1-2H')
    assert (mol.inchi.key_layer_content(C2H2F2_ICH_STEREO_UNKNOWN, 'b')
            == '2-1?')
    assert mol.inchi.key_layer_content(C2H2F2_ICH_NO_STEREO, 'c') == '3-1-2-4'
    assert mol.inchi.key_layer_content(C2H2F2_ICH_NO_STEREO, 'h') == '1-2H'
    assert mol.inchi.key_layer_content(C2H2F2_ICH_NO_STEREO, 'b') is None


def test__inchi__core_parent():
    """ inchi.core_parent
    """
    assert mol.inchi.core_parent(AR_ICH) == AR_ICH
    assert (mol.inchi.core_parent(C2H2F2_ICH)
            == 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H')
    assert (mol.inchi.core_parent(C2H2F2_ICH_NO_STEREO)
            == 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H')
    assert (mol.inchi.core_parent(C2H2F2_ICH_STEREO_UNKNOWN)
            == 'InChI=1/C2H2F2/c3-1-2-4/h1-2H')
    assert (mol.inchi.core_parent(C8H13O_ICH_NO_STEREO)
            == 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3')
    assert (mol.inchi.core_parent(C8H13O_ICH)
            == 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3')


def test__inchi__atom_stereo_elements():
    """ inchi.atom_stereo_elements
    """
    assert mol.inchi.atom_stereo_elements(C8H13O_ICH_NO_STEREO) == ()
    assert mol.inchi.atom_stereo_elements(C8H13O_ICH) == (('8', '-'),)


def test__inchi__bond_stereo_elements():
    """ inchi.bond_stereo_elements
    """
    assert mol.inchi.bond_stereo_elements(C8H13O_ICH_NO_STEREO) == ()
    assert (mol.inchi.bond_stereo_elements(C8H13O_ICH)
            == (('5-3', '-'), ('6-4', '-')))


def test__inchi__has_unknown_stereo_elements():
    """ inchi.has_unknown_stereo_elements
    """
    assert (mol.inchi.has_unknown_stereo_elements(C8H13O_ICH)
            is False)
    assert (mol.inchi.has_unknown_stereo_elements(C8H13O_ICH_PARTIAL_STEREO)
            is True)
    assert (mol.inchi.has_unknown_stereo_elements(C8H13O_ICH_NO_STEREO)
            is True)
    assert (mol.inchi.has_unknown_stereo_elements(C8H13O_ICH_NO_ENANTIOMER)
            is False)


def test__inchi__compatible_stereoisomers():
    """ inchi.compatible_stereoisomers
    """
    assert (mol.inchi.compatible_stereoisomers(C8H13O_ICH_NO_STEREO)
            == C8H13O_ICHS)
    assert mol.inchi.compatible_stereoisomers(C8H13O_ICH) == (C8H13O_ICH,)


def test__inchi__inchi_key():
    """ inchi.inchi_key
    """
    assert mol.inchi.inchi_key(C2H2F2_ICH) == 'WFLOTYSKFUPZQB-OWOJBTEDSA-N'
    assert (mol.inchi.inchi_key(C2H2F2_ICH_NO_STEREO)
            == 'WFLOTYSKFUPZQB-UHFFFAOYSA-N')


def test__inchi__key__first_hash():
    """ inchi.key.first_hash()
    """
    assert (mol.inchi.key.first_hash(
        mol.inchi.inchi_key(C2H2F2_ICH)) == 'WFLOTYSKFUPZQB')
    assert (mol.inchi.key.first_hash(
        mol.inchi.inchi_key(C2H2F2_ICH_NO_STEREO)) == 'WFLOTYSKFUPZQB')
    assert (mol.inchi.key.first_hash(
        mol.inchi.inchi_key(C2H2F2_ICH_STEREO_UNKNOWN)) == 'WFLOTYSKFUPZQB')


def test__inchi__key__second_hash():
    """ inchi.key.second_hash()
    """
    assert (mol.inchi.key.second_hash(
        mol.inchi.inchi_key(C2H2F2_ICH)) == 'OWOJBTED')
    assert (mol.inchi.key.second_hash(
        mol.inchi.inchi_key(C2H2F2_ICH_NO_STEREO)) == 'UHFFFAOY')
    assert (mol.inchi.key.second_hash(
        mol.inchi.inchi_key(C2H2F2_ICH_STEREO_UNKNOWN)) == 'HXYFBOIP')


def test__inchi__key__version_indicator():
    """ inchi.key.version_indicator()
    """
    assert (mol.inchi.key.version_indicator(
        mol.inchi.inchi_key(C2H2F2_ICH)) == 'SA')
    assert (mol.inchi.key.version_indicator(
        mol.inchi.inchi_key(C2H2F2_ICH_NO_STEREO)) == 'SA')
    assert (mol.inchi.key.version_indicator(
        mol.inchi.inchi_key(C2H2F2_ICH_STEREO_UNKNOWN)) == 'NA')


def test__inchi__key__protonation_indicator():
    """ inchi.key.protonation_indicator()
    """
    ich1 = 'InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)'
    ich2 = 'InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)/p-1'
    ich3 = 'InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)/p+1'
    assert (mol.inchi.key.protonation_indicator(mol.inchi.inchi_key(ich1))
            == 'N')
    assert (mol.inchi.key.protonation_indicator(mol.inchi.inchi_key(ich2))
            == 'M')
    assert (mol.inchi.key.protonation_indicator(mol.inchi.inchi_key(ich3))
            == 'O')


def test__inchi__geometry():
    """ inchi.geometry
    """
    # make sure these run
    mol.inchi.geometry(C2H2F2_ICH)
    for ich in C8H13O_ICHS:
        mol.inchi.geometry(ich)
    for ich in RDKIT_FAIL_ICHS:
        mol.inchi.geometry(ich)
    for ich in PYBEL_FAIL_ICHS:
        mol.inchi.geometry(ich)


def test__inchi__connectivity_graph():
    """ inchi.connectivity_graph
    """
    assert (mol.inchi.connectivity_graph(C2H2F2_ICH)
            == ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('F', 0, None),
                 3: ('F', 0, None), 4: ('H', 0, None), 5: ('H', 0, None)},
                {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
                 frozenset({1, 3}): (1, None), frozenset({0, 4}): (1, None),
                 frozenset({1, 5}): (1, None)}))

    # make sure this runs -- the function is self-testing
    for ich in HEPTANE_ICHS:
        mol.inchi.connectivity_graph(ich)


def test__inchi__stereo_graph():
    """ inchi.stereo_graph
    """
    assert (mol.inchi.stereo_graph(C8H13O_ICH)
            == ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
                 3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 0, None),
                 6: ('C', 0, None), 7: ('C', 0, False), 8: ('O', 0, None),
                 9: ('H', 0, None), 10: ('H', 0, None), 11: ('H', 0, None),
                 12: ('H', 0, None), 13: ('H', 0, None), 14: ('H', 0, None),
                 15: ('H', 0, None), 16: ('H', 0, None), 17: ('H', 0, None),
                 18: ('H', 0, None), 19: ('H', 0, None), 20: ('H', 0, None),
                 21: ('H', 0, None)},
                {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
                 frozenset({2, 4}): (1, False), frozenset({3, 5}): (1, False),
                 frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
                 frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
                 frozenset({0, 9}): (1, None), frozenset({0, 10}): (1, None),
                 frozenset({0, 11}): (1, None), frozenset({1, 12}): (1, None),
                 frozenset({1, 13}): (1, None), frozenset({1, 14}): (1, None),
                 frozenset({2, 15}): (1, None), frozenset({16, 3}): (1, None),
                 frozenset({17, 4}): (1, None), frozenset({18, 5}): (1, None),
                 frozenset({19, 6}): (1, None), frozenset({20, 6}): (1, None),
                 frozenset({21, 7}): (1, None)}))

    # make sure this runs -- the function is self-testing
    ichs = numpy.random.choice(HEPTANE_ICHS, 10)
    for ich in ichs:
        try:
            mol.inchi.stereo_graph(ich)
        except NotImplementedError:
            pass


if __name__ == '__main__':
    # test__smiles__inchi()
    # test__inchi__smiles()
    # test__inchi__recalculate()
    # test__inchi__is_closed()
    # test__inchi__prefix()
    # test__inchi__version()
    # test__inchi__formula_layer()
    # test__inchi__key_layer()
    # test__inchi__key_layer_content()
    # test__inchi__core_parent()
    # test__inchi__atom_stereo_elements()
    # test__inchi__bond_stereo_elements()
    # test__inchi__has_unknown_stereo_elements()
    # test__inchi__compatible_stereoisomers()
    # test__inchi__key__first_hash()
    # test__inchi__key__second_hash()
    test__inchi__key__version_indicator()
    test__inchi__key__protonation_indicator()
    # test__inchi__geometry()
    # test__inchi__connectivity_graph()
    # test__inchi__stereo_graph()
