""" graph geometry library

For generating heuristic coordinates and z-matrices from graphs
"""
import numpy
import scipy.optimize
from qcelemental import constants as qcc
from qcelemental import periodictable as pt
import automol.zmat
from automol.graph._res import resonance_dominant_atom_hybridizations
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import atom_symbols
from automol.graph._graph import longest_chain
from automol.graph._graph import rings_sorted_atom_keys


# bond distances
XY_DIST = 1.5       # angstroms
XH_DIST = 1.1       # angstroms

# bond angles
TET_ANG = 109.4712  # degrees
TRI_ANG = 120.      # degrees
LIN_ANG = 180.      # degrees

# dihedral angles
CIS_DIH = 0.        # degrees
TRA_DIH = 180.      # degrees

# miscellaneous
RIT_ANG = 90.       # degrees


def heuristic_bond_distance(gra, key1, key2, check=True):
    """ heuristic bond distance
    """
    if check:
        assert key1 in atom_neighbor_keys(gra)[key2]

    sym_dct = atom_symbols(gra)
    sym1 = sym_dct[key1]
    sym2 = sym_dct[key2]

    if pt.to_Z(sym1) == 1 or pt.to_Z(sym2) == 1:
        dist = XH_DIST
    else:
        dist = XY_DIST

    return dist


def heuristic_bond_angle(gra, key1, key2, key3, check=True):
    """ heuristic bond angle
    """
    if check:
        assert {key1, key3} <= set(atom_neighbor_keys(gra)[key2])

    hyb_dct = resonance_dominant_atom_hybridizations(gra)
    hyb2 = hyb_dct[key2]
    if hyb2 == 3:
        ang = TET_ANG
    elif hyb2 == 2:
        ang = TRI_ANG
    else:
        assert hyb2 == 1
        ang = LIN_ANG

    return ang


def chain_zmatrix(gra, chain_keys):
    """ build a chain of heavy atoms
    """
    chain_iter = iter(chain_keys)

    sym_dct = atom_symbols(gra)     # graph keys -> atomic symbols
    row_dct = {}                    # graph keys -> z-matrix rows
    zma = ()                        # empty z-matrix

    key3 = row3 = r34 = None
    key2 = row2 = a234 = None
    key1 = row1 = d1234 = None

    for key4 in chain_iter:
        row_dct[key4] = automol.zmat.count(zma)

        if key3 is not None:
            row3 = row_dct[key3]
            r34 = heuristic_bond_distance(gra, key3, key4)

        if key2 is not None:
            row2 = row_dct[key2]
            a234 = heuristic_bond_angle(gra, key2, key3, key4)

        if key1 is not None:
            row1 = row_dct[key1]
            d1234 = TRA_DIH

        zma = automol.zmat.add_atom(zma, sym_dct[key4],
                                    key_row=[row3, row2, row1],
                                    val_row=[r34, a234, d1234])

        # now, shift the keys for the next one up
        key1, key2, key3 = key2, key3, key4

    chain_rows = tuple(map(row_dct.__getitem__, chain_keys))
    return zma, chain_rows


def ring_arc_bond_angle(num, end_dist=XY_DIST, bond_dist=XY_DIST):
    """ find the angle subtended by a ring arc (in degrees)

    We will form an arc of n atoms such that the ends are a distance r_e
    apart and each pair of neighboring atoms is a distance r_b apart.
    Let the total angle subtended by the arc be theta and let the circle
    radius be R.
    Furthermore, let the arc angle subtended between two neighboring atoms be
    alpha, so that
        alpha = theta/(n-1)
    Let us consider two isosceles triangles: the one formed by two neighboring
    atoms on the arc with the circle center (triangle N), and the one formed by
    the arc ends with the circle center (triangle E).

    By bisecting triangle N, we find that
        sin(alpha/2) = (r_b/2)/R = r_b/(2*R)
    By bisecting triangle E, we find that
        +-sin(theta/2) = (r_e/2)/R = r_e/(2*R)
    where the sign flips when theta > pi.
    Dividing these two equations gives us.
        sin(alpha/2)/sin(theta/2) = r_b/r_e
    By substituting in the expression for alpha in terms of theta, we can solve
    for theta by finding the first positive root of the following equation.
        sin(theta/(2*(n-1)))/sin(theta/2) - r_b/r_e = 0
    Note that, in solving this equation numerically, we must avoid the
    singularity at theta = 0.

    Once we know theta, we can find the bond angle for atoms along the arc as
    follows.
    Let beta be the base angle of the triangle N, and note that alpha is the
    vertex angle.
    Geometrically, then, beta is the bisection of the bond angle, and therefore
    the bond angle is
        bond_angle = 2*beta.
    Furthermore, by the angle sum formula for triangle N, we have that
        alpha + 2*beta = 180 degrees (pi radians)
    It therefore follows that
        bond_angle = 180 - alpha = 180 - theta/(n-1)
    """

    def _f(theta):
        if numpy.abs(theta) < 0.001:
            theta = 0.001

        rhs = numpy.sin(theta/(2*(num-1)))/numpy.sin(theta/2)
        lhs = bond_dist / end_dist
        return rhs - lhs

    res_obj = scipy.optimize.root_scalar(_f, method='brentq',
                                         bracket=[0.01, 2*numpy.pi])

    assert res_obj.converged
    # convert from radians to degrees
    theta = res_obj.root * qcc.conversion_factor('radian', 'degree')

    return 180. - theta/(num-1)


def ring_zmatrix(gra, ring_keys, bond_dist=XY_DIST, end_dist=XY_DIST):
    """ build an arc (or ring) of heavy atoms

    :param bond_dist: bond distances between neighboring atoms in the ring
    :type bond_dist: float
    :param end_dist: if this is an arc connecting two points of a larger ring,
        this sets the distance between the arc ends
    :type end_dist: float
    """

    num = len(ring_keys)
    bond_ang = ring_arc_bond_angle(num, end_dist=end_dist, bond_dist=bond_dist)
    print(bond_ang)

    print(ring_keys)


if __name__ == '__main__':
    import automol

    # chain
    ICH = automol.smiles.inchi('CCCCCC')
    GRA = automol.inchi.graph(ICH)
    # print(automol.graph.string(GRA, one_indexed=False))
    CHAIN_KEYS = longest_chain(GRA)
    print(CHAIN_KEYS)
    # ZMA, CHAIN_ROWS = chain_zmatrix(GRA, CHAIN_KEYS)

    # print(CHAIN_ROWS)
    # print(automol.zmat.string(ZMA))
    # GEO = automol.zmat.geometry(ZMA)
    # print(automol.geom.string(GEO))

    # ring
    ICH = automol.smiles.inchi('C12CC(C2)CC1')
    GRA = automol.inchi.graph(ICH)
    print(automol.graph.string(GRA, one_indexed=False))
    RING_KEYS_LST = sorted(rings_sorted_atom_keys(GRA), key=len)
    ring_zmatrix(GRA, RING_KEYS_LST[1])
