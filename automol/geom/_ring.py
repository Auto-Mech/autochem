""" get lvl 4
"""

import numpy as np
from phydat import phycon

from automol.geom._conv import (
    graph,
)
from automol.geom.base import central_angle, subgeom
from automol.graph import base as graph_base

ATHRESH = 80.0 * phycon.DEG2RAD


def all_rings_angles_reasonable(geo, rings_atoms, thresh=ATHRESH):
    """Assess if all angles of ring are cool"""

    condition = True
    for ring_atoms in rings_atoms:
        condition = ring_angles_reasonable(geo, ring_atoms, thresh=thresh)

    return condition


def ring_angles_reasonable(geo, ring_atoms, thresh=ATHRESH):
    """Assess whether any of the angles of the prospective rings of a geometry
    are bent at unphysical angles.

    :param geo: molecular geometry
    :type geo: automol.geom object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    :param thresh: angular threshold for large angles
    :type thresh: float
    :rtype: bool
    """

    condition = True
    # Lower threshold for 4 and 5 membered rings if not using relaxed thresh already
    if len(ring_atoms) < 6 and thresh > 70.: thresh *= 0.8
    for i, ring_atom in enumerate(ring_atoms):
        _atoms = [ring_atom, ring_atoms[i - 1], ring_atoms[i - 2]]
        cangle = central_angle(geo, *_atoms, degree=False)
        if cangle < thresh:
            condition = False
            break

    return condition


def ring_fragments_geometry(geo, rings_atoms=None, ngbs=None):
    """Fragment out the ring and its neighbors in a geometry.(?)

    :param geo: molecular geometry
    :type geo: automol.geom object
    """
    #adl TODO call with fix_hyper = False to prevent "assertion" error when calling fram_samp_geo
    gra = graph(geo)  
    if rings_atoms is None:
        rings_atoms = graph_base.rings_atom_keys(gra)
    if ngbs is None:
        ngbs = graph_base.atoms_sorted_neighbor_atom_keys(gra)

    ring_idxs = []
    ret = None
    for ring_atoms in rings_atoms:
        for ring_atom in ring_atoms:
            ring_ngbs = ngbs[ring_atom]
            if ring_atom not in ring_idxs:
                ring_idxs.append(ring_atom)
            for ngb in ring_ngbs:
                if ngb not in ring_idxs:
                    ring_idxs.append(ngb)
    if ring_idxs:
        ret = subgeom(geo, ring_idxs)

    return ret


# adl added functions for Cremer Pople parameters
# Inspired from Chan et al. 2021 https://doi.org/10.1021/acs.jcim.0c01144 
def translate_to_ring_center(coords):
    """Translate reference system to center of the ring

    :param coords: cartesian coordinates of ring atoms
    :type coords: list of lists [[x1,y1,z1],[x2,y2,z2],...]

    :output new_coords: translated cartesian coordinates of ring atoms
    :type new_coords: np.array 
    """
    new_coords = coords - coords.mean(axis=0)
    return new_coords

def mean_ring_plane(coords):
    """Compute the mean ring plane

    :param coords: cartesian coordinates of ring atoms, origin on ring center
    :type coords: np.array

    :output R1,R2: 3d vectors defining the ring mean plane
    :type R1,R2: np.array  
    """
    n_ring_atoms = coords.shape[0] # ring size
    R1 = np.dot(np.sin(2*np.pi*np.arange(0,n_ring_atoms)/n_ring_atoms),coords)
    R2 = np.dot(np.cos(2*np.pi*np.arange(0,n_ring_atoms)/n_ring_atoms),coords)
    return R1, R2

def normal_to_ring_plane(coords):
    """Compute the normalized normla to the mean ring plane

    :param coords: cartesian coordinates of ring atoms, origin on ring center
    :type coords: np.array

    :output n: normalized vector perpendicular to R1,R2 plane
    :type n: np.array  
    """
    R1,R2 = mean_ring_plane(coords)
    n = np.cross(R1,R2)/np.linalg.norm(np.cross(R1,R2))
    return n

def get_displacement(coords):
    """Compute atoms displacement along n direction with respect to ring plane

    :param coords: cartesian coordinates of ring atoms, origin on ring center
    :type coords: np.array

    :output z: vector defining the atoms dispalcements with respect to the ring mean plane
    :type z: np.array  
    """
    n = normal_to_ring_plane(coords)
    z = np.dot(coords,n)
    return z

def cremer_pople_params(coords):
    """Computes Cremer-Pople parameters for conformational analysis of ring structures.

    :param coords: cartesian coordinates of ring atoms
    :type coords: list of lists [[x1,y1,z1],[x2,y2,z2],...]

    :output amplitude,angle: set of puckering parameters
    :type amplitude,angle: tuple of lists ([q2,q3,...],[phi2,phi3...])
    """
    # cremer_pople_params function:
    coords_array = np.array(coords)
    coords_array = translate_to_ring_center(coords_array)

    n_atoms = coords_array.shape[0]  # number of atoms in the ring
    z = get_displacement(coords_array)
    if n_atoms%2 == 0: # n_atoms even
        m = range(2,int((n_atoms/2)))
        cos_component = [np.dot(z,np.cos(2*np.pi*k*np.arange(0,n_atoms)/n_atoms)) for k in m]
        sin_component = [np.dot(z,np.sin(2*np.pi*k*np.arange(0,n_atoms)/n_atoms)) for k in m]
        qcos = np.sqrt(2/n_atoms)*np.array(cos_component)
        qsin = -np.sqrt(2/n_atoms)*np.array(sin_component)
        q = np.sqrt(qsin**2 + qcos**2)
        amplitude = np.append(q, (1/np.sqrt(n_atoms))*np.dot(z,np.cos(np.arange(0,n_atoms)*np.pi)).sum()).tolist()
        angle = np.arctan2(qsin,qcos).tolist()
    else: # n_atoms odd
        m = range(2,int((n_atoms-1)/2)+1)
        cos_component = [np.dot(z,np.cos(2*np.pi*k*np.arange(0,n_atoms)/n_atoms)) for k in m]
        sin_component = [np.dot(z,np.sin(2*np.pi*k*np.arange(0,n_atoms)/n_atoms)) for k in m]
        qcos = np.sqrt(2/n_atoms)*np.array(cos_component)
        qsin = -np.sqrt(2/n_atoms)*np.array(sin_component)
        amplitude = np.sqrt(qsin**2 + qcos**2).tolist()
        angle = np.arctan2(qsin,qcos).tolist()

    return (amplitude, angle), z.tolist()
    





