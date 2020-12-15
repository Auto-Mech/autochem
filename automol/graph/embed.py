""" geometry embedding using the distance geometry algorithm

Blaney, J. M.; Dixon, J. S. “Distance Geometry in Molecular Modeling”.
Reviews in Computational Chemistry; VCH: New York, 1994.
"""
import itertools
import more_itertools as mit
import numpy
import qcelemental as qce
from qcelemental import periodictable as pt
from automol.graph._graph_base import string
from automol.graph._graph import atom_symbols
from automol.graph._graph import atom_keys
from automol.graph._graph import atom_shortest_paths
from automol.graph._graph import atom_neighbor_keys


# bond distances
XY_DIST = 1.5       # angstroms
XH_DIST = 1.1       # angstroms


def heuristic_bond_distance(gra, key1, key2, check=True):
    """ heuristic bond distance (in angstroms)
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


def van_der_waals_radius(gra, key):
    """ van der waals radius for an atom in the graph (in angstroms)
    """
    sym_dct = atom_symbols(gra)
    sym = sym_dct[key]
    rad = qce.vdwradii.get(sym, units='angstrom')
    return rad


def closest_approach(gra, key1, key2):
    """ closest approach between atoms, based on their van der Waals radii

    Warning: The scaling factor on the van der waals radii was arbitrarily
    chosen based on limited tests and may need to be lowered
    """
    vdw_scaling_factor = 0.9
    dist = (van_der_waals_radius(gra, key1) +
            van_der_waals_radius(gra, key2)) * vdw_scaling_factor
    return dist


def upper_distance_bound(gra, path):
    """ upper distance bound between two ends of a path

    :param path: the shortest path between two atoms
    :type path: list or tuple
    """
    # if the path is 0, the atoms are disconnected and could be arbitrarily far
    # apart
    if len(path) == 1:
        dist = 0
    elif len(path) == 2:
        dist = heuristic_bond_distance(gra, *path)
    # otherwise, just do the sum of the distances between atoms along the path
    else:
        assert len(path) > 2
        print('HERE')
        print(path)
        dist = sum((heuristic_bond_distance(gra, *bond)
                    for bond in mit.windowed(path, 2)))

    return dist


def lower_distance_bound(gra, path):
    """ lower distance bound between two ends of a path

    :param path: the shortest path between two atoms
    :type path: list or tuple
    """
    if len(path) == 1:
        dist = 0
    elif len(path) == 2:
        dist = heuristic_bond_distance(gra, *path)
    # otherwise, just use the VDW radii for the two atoms
    else:
        assert len(path) > 2
        dist = closest_approach(gra, path[0], path[-1])

    return dist


def distance_bounds_matrix(gra, keys):
    """ initial distance bounds matrix
    """
    assert set(keys) == set(atom_keys(gra))

    sp_dct = atom_shortest_paths(gra)

    natms = len(keys)
    bmat = numpy.zeros((natms, natms))
    for i, j in itertools.combinations(range(natms), 2):
        if j in sp_dct[i]:
            path = sp_dct[i][j]
            bmat[i, j] = upper_distance_bound(gra, path)
            bmat[j, i] = lower_distance_bound(gra, path)
        else:
            # they are disconnected
            bmat[i, j] = 999
            bmat[j, i] = closest_approach(gra, keys[i], keys[j])

        assert bmat[j, i] <= bmat[i, j], (
            "Lower bound exceeds upper bound. Something needs to be fixed!\n"
            "{}\npath: {}\n"
            .format(string(gra, one_indexed=False), str(path)))

    return bmat


if __name__ == '__main__':
    import automol
    ICH = automol.smiles.inchi('C1CC1.C')
    GRA = automol.inchi.graph(ICH)
    GRA = automol.graph.explicit(GRA)
    KEYS = sorted(automol.graph.atom_keys(GRA))
    BMAT = distance_bounds_matrix(GRA, KEYS)
    print(BMAT)
