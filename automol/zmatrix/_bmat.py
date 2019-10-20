"""
IN PROGRESS:Transliterating Carlo's routine from Fortran

Form the B-Matrix and C-Matrix used to convert the coordinates

Calcualtes all of the derivaties via finite-difference

define starting xyz geometry.
convention: atom 1 is is 0 0 0
atom 2 bd 0 0
atom 3 on xy plane
"""

import numpy as np

NATOMS = 10  # maybe need, don't know
INT_COORDS = ''
DELTAX = 0.01
DELTAY = 0.01
DELTAZ = 0.01






def compute_bmat(natoms, coords, deltax, deltay, deltaz):
    """ compute the bmatrix by central difference
        where B_ik = dq_i / dx_k
    """

    b_mat = np.zeros(3*natoms, 3*natoms)
    for j in range(3):
        for k in range(3):

            # perturb x + dx and x - dx
            xpert_xp = 1
            xpert_xn = 1
            _perturb_coordinates(coords, jpert, delta)

            # perturb y + dy and y - dy
            xpert_yp = 1
            xpert_yn = 1
            _perturb_coordinates(coords, jpert, delta)

            # perturb z + dz and z - dz
            xpert_zp = 1
            xpert_zn = 1
            _perturb_coordinates(coords, jpert, delta)

        # Now calculate the jk component C-Matrix
        _calculate_bmat_k_component(b_mat, coords, j, j*k,
                                    x_pert_pp, x_pert_pn,
                                    x_pert_np, x_pert_nn)

        # now update iangsub1 bmat component (whatever this is)
        b_mat = _update_bmat(bmat, coords)

    return b_mat


def compute_cmat(natoms, coords, deltax, deltay, deltaz):
    """ compute the bmatrix by central difference
        where C_ijk = d2q_i / (dx_j.dx_k)
    """

    c_mat = np.zeros(3*natoms, 3*natoms, 3*natoms)
    for j in range(3):
        for k in range(3):
            # perturb xj + dxj and xk + dxk
            x_pert_pp = _perturb_coordinates(coords, jpert, kpert, d1, d2)

            # perturb xj - dxj and yk + dyk
            x_pert_np = _perturb_coordinates(coords, jpert, kpert, d1, d2)

            # perturb xj + dxj and yk - dyk
            x_pert_pn = _perturb_coordinates(coords, jpert, kpert, d1, d2)

            # perturb xj - dxj and xk - dxk
            x_pert_nn = _perturb_coordinates(coords, jpert, kpert, d1, d2)

            # Now calculate the jk component C-Matrix
            _calculate_cmat_k_component(c_mat, coords, j, j*k,
                                        x_pert_pp, x_pert_pn,
                                        x_pert_np, x_pert_nn)

    return c_mat


def _perturb_coordinates(coords, jpert, delta1, kpert=None, delta2=None):
    """ Generate coordinates that have been perturbed
    """
    coords[jpert] += delta1
    coords[kpert] += delta2
    #           call update_zmat(natom,natomt,intcoor,bislab,ibconn,
    # $     iaconn,idconn,bname,anname,dname,atname,cooxpp,cooypp,
    # $     coozpp,xintpp,tauopt,ntau,idummy,ilin_fr,aconnt,bconnt,
    # $     dconnt,atomlabel,ifilu)

    return coords


def _calculate_bmat_k_component(b_mat, j_idx, coords, delta,
                                x_pert_p, x_pert_n):
    """ Calculate one nine components of B_ij for given __
    """

    for i, coord in enumerate(coords):
        if abs(xpert_p[i] - xpert_np[i]) > 300.0:
            if xpert_n[i] < 0.0:
                xpert_n[i] += 360.0
            elif xpert_n[i] > 0.0:
                xpert_n[i] -= 360.0
            if abs(xpert_p[i] - xpert_n[i]) > 300.0:
                raise ValueError(
                    'something did not work here: k, j coord', kind, jind, i)
        b_mat[i, j_idx] = (
            ((xpert_p[i] - xpert_n[i]) / 2.0) * (1.0 / delta)
        )

    return b_mat


def _calculate_cmat_k_component(c_mat, k_idx, coords, delta1, delta2,
                                x_pert_pp, x_pert_pn, x_pert_np, x_pert_nn):
    """ Calculate one nine components of C_ijk for given j
    """

    for i, coord in enumerate(coords):

        if abs(xpert_pp[i] - xpert_np[i]) > 300.0:
            if xpert_pp[i] < 0.0:
                xpert_pp[i] += 360.0
            elif xpert_pp[i] > 0.0:
                xpert_pp[i] -= 360.0
            if abs(xpert_pp[i] - xpert_np[i]) > 300.0:
                raise ValueError(
                    'something did not work here: k, j coord',
                    kind, jind, i)

        if abs(xpert_np[i] - xpert_np[i]) > 300.0:
            if xpert_pn[i] < 0.0:
                xpert_pn[i] += 360.0
            elif xpert_pn[i] > 0.0:
                xpert_pn[i] -= 360.0
            if abs(xpert_pp[i] - xpert_pn[i]) > 300.0:
                raise ValueError(
                    'something did not work here: k, j coord',
                    kind, jind, i)

        if abs(xpert_np[i] - xpert_nn[i]) > 300.0:
            if xpert_nn[i] < 0.0:
                xpert_nn[i] += 360.0
            elif xpert_nn[i] > 0.0:
                xpert_nn[i] -= 360.0
                if abs(xpert_np[i] - xpert_nn[i]) > 300.0:
                    raise ValueError(
                        'something did not work here: k, j coord',
                        kind, jind, i)

        c_mat[i, j_idx, k_idx] = (
            xpert_pp[i] - xpert_np[i] - xpert_pn[i] +
            (xpert_nn[i] / 4.0) * (1.0 / deltax) * (1.0 / deltaz)
        )

    return c_mat


if __name__ == '__main__':
    b_mat = compute_bmat(NATOMS, COORDS, DELTAX, DELTAY, DELTAZ)
    c_mat = compute_cmat(NATOMS, COORDS, DELTAX, DELTAY, DELTAZ)
