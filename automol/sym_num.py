""" experimentalsymmetry number codes
"""
import itertools
import numpy
import sympy.combinatorics as spc
from automol import util
import automol.geom
from automol.convert._pyx2z import to_oriented_geometry


# symmetry number codes
def external_symmetry_factor1(geo):
    """ obtain external symmetry factor for a geometry using x2z
    """
    # Get initial external symmetry number
    if automol.geom.is_atom(geo):
        ext_sym_fac = 1.
    else:
        oriented_geom = to_oriented_geometry(geo)
        ext_sym_fac = oriented_geom.sym_num()
    return ext_sym_fac


def external_symmetry_factor2(geo, thresh=1e-3):
    """ calculate the external symmetry factor
    """
    sym_num = 1

    if automol.geom.is_linear(geo):
        print("LINEAR")
    elif not automol.geom.is_atom(geo):
        idxs = range(automol.geom.count(geo))
        orig = (0., 0., 0.)

        # First, find two atoms that aren't at the origin and aren't on a line
        # through the origin
        ref_geo = automol.geom.mass_centered(geo)
        for idx_pair in itertools.permutations(idxs, r=2):
            xyzs = automol.geom.coordinates(ref_geo, idxs=idx_pair)
            norms = list(map(numpy.linalg.norm, xyzs))
            if all(norm > thresh for norm in norms):
                ang = util.vec.central_angle(orig, *xyzs)
                if not (numpy.abs(ang) < thresh or
                        numpy.abs(ang - numpy.pi) < thresh):
                    ref_idx_pair = idx_pair
                    break

        # This pair of atoms, together with the center of mass, will serve as
        # our "reference points" that we can align with other orientations as
        # we identify symmetries
        ref_xyzs = automol.geom.coordinates(ref_geo, idxs=ref_idx_pair)
        ref_syms = automol.geom.symbols(ref_geo, idxs=ref_idx_pair)

        ref_xyzs_with_orig = (orig,) + ref_xyzs
        ref_dist_mat = util.vec.distance_matrix(ref_xyzs_with_orig)

        # Now, loop over pairs. If they have the same symbols and the same
        # distance matrix with the origin as the reference pair, rotate the
        # geometry about the origin to superimpose them onto the reference
        # pair. If the rotated geometry is exactly aligned with the reference
        # geometry, then we've identified a new symmetry to add to our list.
        perms = []
        for idx_pair in itertools.permutations(idxs, r=2):
            geo = ref_geo
            syms = automol.geom.symbols(geo, idxs=idx_pair)
            if idx_pair != ref_idx_pair and syms == ref_syms:
                xyzs = automol.geom.coordinates(geo, idxs=idx_pair)
                xyzs_with_orig = (orig,) + xyzs
                dist_mat = util.vec.distance_matrix(xyzs_with_orig)
                if numpy.allclose(dist_mat, ref_dist_mat):
                    rot_mat = util.mat.superimposition(xyzs_with_orig,
                                                       ref_xyzs_with_orig)
                    geo = automol.geom.transform_by_matrix(geo, rot_mat)
                    perm = automol.geom.permutation(
                        geo, ref_geo, thresh=thresh)
                    if perm is not None:
                        perms.append(perm)

        sym_num = _symmetry_number(perms)

    sym_fac = sym_num
    return sym_fac


def external_symmetry_factor3(geo, thresh=1e-3):
    """ calculate the external symmetry factor
    """
    sym_num = 1

    if automol.geom.is_linear(geo):
        print("LINEAR")
    elif not automol.geom.is_atom(geo):
        idxs = range(automol.geom.count(geo))

        # First, find three atoms that aren't on a line
        ref_geo = automol.geom.mass_centered(geo)
        for idx_trip in itertools.permutations(idxs, r=3):
            if not automol.geom.is_linear(
                    automol.geom.from_subset(geo, idxs=idx_trip)):
                ref_idx_trip = idx_trip
                break

        # These three atoms serve as reference points for aligning orientations
        # as we identify symmetries
        ref_xyzs = automol.geom.coordinates(ref_geo, idxs=ref_idx_trip)
        ref_syms = automol.geom.symbols(ref_geo, idxs=ref_idx_trip)

        ref_dist_mat = util.vec.distance_matrix(ref_xyzs)

        # Now, loop over triplets in the geometry. If they have the same
        # symbols and distance matrix as the refernce triplet, rotate the
        # geometry to superimpose them onto the reference triplet. If the
        # rotated geometry is a permutation of the present geometry, add the
        # permutation to our list of symmetries.
        perms = []
        for idx_trip in itertools.permutations(idxs, r=3):
            geo = ref_geo
            syms = automol.geom.symbols(geo, idxs=idx_trip)
            if syms == ref_syms:
                xyzs = automol.geom.coordinates(geo, idxs=idx_trip)
                dist_mat = util.vec.distance_matrix(xyzs)
                if numpy.allclose(dist_mat, ref_dist_mat):
                    rot_mat = util.mat.superimposition(xyzs, ref_xyzs)
                    geo = automol.geom.transform_by_matrix(geo, rot_mat)
                    perm = automol.geom.permutation(
                        geo, ref_geo, thresh=thresh)
                    if perm is not None:
                        perms.append(perm)

        sym_num = _symmetry_number(perms)

    sym_fac = sym_num
    return sym_fac


def _symmetry_number(perms):
    """ calculate the symmetry number from a list of permutations
    """
    perm_objs = list(map(spc.permutations.Permutation, perms))
    perm_group = spc.perm_groups.PermutationGroup(perm_objs)
    sym_num = perm_group.order()
    return sym_num


if __name__ == '__main__':
    # GEO = (('C', (1.43210415746, -0.7681543652, 0.212793208918)),
    #        ('C', (-0.0393296328, 1.6305747247, -0.161072379111)),
    #        ('C', (-1.3927747481, -0.8624201904, -0.051721285230)),)

    GEO = (('C', (0.0, 0.0, 0.0)),
           ('H', (0.0, 0.0, 2.0786987380036113)),
           ('H', (0.0, 1.9598159649150293, -0.6928995793345372)),
           ('H', (1.69725041235873, -0.97990798245751, -0.69289957933454)),
           ('H', (-1.69725041235873, -0.97990798245751, -0.69289957933454)))
    print(external_symmetry_factor2(GEO))
