""" get lvl 4
"""

import os
import subprocess
import numpy as np

from phydat import phycon

from ..graph import base as graph_base
from ._1conv import (
    graph,
)
from .base import (
    central_angle,
    dihedral_angle,
    subgeom,
    from_xyz_trajectory_string,
    xyz_string,
    string,
)

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
    # Lower threshold for 4 and 5 membered rings
    # if not using relaxed thresh already
    if len(ring_atoms) < 6 and thresh > 70.:
        thresh *= 0.8
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

def ring_only_geometry(geo, rings_atoms=None):
    """Fragment out the rings in a geometry.
    How to deal with geos in which rings are not connected??

    :param geo: molecular geometry
    :type geo: automol.geom object
    """

    if rings_atoms is None:
        gra = graph(geo)
        rings_atoms = graph_base.rings_atom_keys(gra)

    ring_idxs = []
    ret = None
    for ring_atoms in rings_atoms:
        for ring_atom in ring_atoms:
            if ring_atom not in ring_idxs:
                ring_idxs.append(ring_atom)

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
    :output z: vector of atoms dispalcements with respect to the ring mean plane
    :type z: np.array  
    """
    n = normal_to_ring_plane(coords)
    z = np.dot(coords,n)
    return z

def cremer_pople_params(coords):
    """Computes Cremer-Pople parameters for conformational analysis of ring structures

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
        cos_component = [np.dot(z,np.cos(
            2*np.pi*k*np.arange(0,n_atoms)/n_atoms)) for k in m]
        sin_component = [np.dot(z,np.sin(
            2*np.pi*k*np.arange(0,n_atoms)/n_atoms)) for k in m]
        qcos = np.sqrt(2/n_atoms)*np.array(cos_component)
        qsin = -np.sqrt(2/n_atoms)*np.array(sin_component)
        q = np.sqrt(qsin**2 + qcos**2)
        amplitude = np.append(q, (1/np.sqrt(n_atoms))*np.dot(
            z,np.cos(np.arange(0,n_atoms)*np.pi)).sum()).tolist()
        angle = np.arctan2(qsin,qcos).tolist()
    else: # n_atoms odd
        m = range(2,int((n_atoms-1)/2)+1)
        cos_component = [np.dot(z,np.cos(
            2*np.pi*k*np.arange(0,n_atoms)/n_atoms)) for k in m]
        sin_component = [np.dot(z,np.sin(
            2*np.pi*k*np.arange(0,n_atoms)/n_atoms)) for k in m]
        qcos = np.sqrt(2/n_atoms)*np.array(cos_component)
        qsin = -np.sqrt(2/n_atoms)*np.array(sin_component)
        amplitude = np.sqrt(qsin**2 + qcos**2).tolist()
        angle = np.arctan2(qsin,qcos).tolist()

    return (amplitude, angle), z.tolist()


def dbscan(geos,rings_atoms, eps, min_samples=1):
    """Density based clustering algorithm
    :param geos: list of geometries
    :type geos: list of automol.geom objects
    :param rings_atoms: list of all atoms in rings
    :type rings_atoms: list of lists of ints
    :output unique_zmas: Z-matrices of unique samples
    :type unique_zmas: list of automol.zmat objects
    :param features: array of features
    :type features: np.array
    :param eps: pradius of a neighborhood centered on a given point
    :type eps: cluster
    :param min_samples: minimum number of points in cluster
    :type min_samples: int
    :output unique_geos: list of unique geometries from cliustering
    :type unique_geos: list of automol.geom objects
    """
    # Possible normalizations, worsened the performance for the puckering
    # def z_score_normalize(features):
    #     """z-score normalization (0 mean, 1 std)

    #     :param features: array of features
    #     :type features: np.array
    #     :output normalized_features: array of normalized features
    #     :type normalized_features: np.array
    #     """
    #     mean = np.mean(features, axis=0)
    #     std = np.std(features, axis=0)
    #     normalized_features = (features - mean) / std
    #     return normalized_features

    # def min_max_normalize(features):
    #     """min-max normalization (all values between 0 and 1)

    #     :param features: array of features
    #     :type features: np.array
    #     :output normalized_features: array of normalized features
    #     :type normalized_features: np.array
    #     """
    #     min_val = np.min(features, axis=0)
    #     max_val = np.max(features, axis=0)
    #     normalized_features = (features - min_val) / (max_val - min_val)
    #     return normalized_features

    def clustering(features, eps, min_samples):

        def find_neighbors(features,point):
            distances = np.linalg.norm(features - features[point], axis=1)
            return np.asarray(distances <= eps).nonzero()[0]

        num_points = features.shape[0]
        labels = np.full(num_points, 0)  # Initialize all labels as 0 (unclassified)
        cluster_id = 1

        for point in range(num_points):
            print(f"Working on point {point}, initial label is {labels[point]}")
            if labels[point] != 0:  # Skip if already classified
                continue

            # Find neighbors
            neighbors = find_neighbors(features,point)
            print(f"Initial neighbors: {neighbors}")
            if len(neighbors) < min_samples:
                labels[point] = -1  # Mark as noise
            else:
                labels[point] = cluster_id
                i = 0
                while i < len(neighbors):
                    neighbor_point = neighbors[i]
                    if labels[neighbor_point] == -1:  # Noise point in cluster
                        labels[neighbor_point] = cluster_id
                    elif labels[neighbor_point] == 0:  # New unvisited point
                        labels[neighbor_point] = cluster_id
                        # Find more neighbors of this point
                        new_neighbors = find_neighbors(features,neighbor_point)
                        if len(new_neighbors) >= min_samples:
                            neighbors = np.concatenate((neighbors, new_neighbors))
                    i += 1
                cluster_id += 1
            print(f"Now label is {labels[point]}")

        return labels

    sub_geos = [ring_only_geometry(geoi,rings_atoms) for geoi in geos]
    subgeo_strings = [xyz_string(geoi, comment="  ") for geoi in sub_geos]

    # Calculate RMSD for each pair of molecules
    import rdkit
    mols = [rdkit.Chem.rdmolfiles.MolFromXYZBlock(geoi) for geoi in subgeo_strings]
    num_mols = len(mols)
    rmsd_matrix = np.zeros((num_mols, num_mols))
    for i in range(num_mols):
        for j in range(i+1, num_mols):
            rmsd = rdkit.Chem.AllChem.GetBestRMS(mols[i], mols[j])
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd

    # Compute displacements from mean ring planes
    z_list = []
    for geoi in geos:
        z_local = []
        for ring_atoms in rings_atoms:
            geo_string = string(geoi, angstrom=True)
            geo_list = [ [float(x) for x in line.split(
                        )[1:]] for line in geo_string.split('\n') ]
            coords = [xyz for i,xyz in enumerate(geo_list) if i in ring_atoms]
            coords = np.array(coords)
            coords = translate_to_ring_center(coords)
            z = get_displacement(coords)
            z_local.extend(z)
        z_list.append(z_local)
    input_z = np.array(z_list)

    # Compute dihedrals of rings atoms
    dih_list = []
    for geoi in geos:
        dih_local = []
        for ring_atoms in rings_atoms:
            for i in range(len(ring_atoms)):
                idxs = [ring_atoms[j%len(ring_atoms)] for j in range(i,i+4)]
                dih = dihedral_angle(geoi,*idxs)
                if dih > np.pi:
                    dih -= 2*np.pi
                elif dih < -np.pi:
                    dih += 2*np.pi
                dih_local.append(dih)
        dih_list.append(dih_local)
    input_dih = np.array(dih_list)

    # Clustering with DBSCAN algorithm
    #input_z = min_max_normalize(input_z)
    features = np.hstack((input_dih,input_z,rmsd_matrix))
    labels = clustering(features, eps=eps, min_samples=min_samples)
    print("labels: ",labels)

    unique_geos = []
    visited_labels = set()
    for label,geoi in zip(labels,geos):
        if label not in visited_labels:
            visited_labels.add(label)
            unique_geos.append(geoi)

    return unique_geos


def checks_with_crest(filename,spc_info):
    """Performs checks with crest on ring geometries to determine unique sampling 
    points for puckering protocol

    :param filename: input file containing geometries, molden style
    :type filename: string
    :param spc_info: input file containing geometries, molden style
    :type spc_info: automech data structure
    :output output unique_geos: list of unique geometries from cliustering
    :type unique_geos: list of automol.geom objects
    """
    # Setup crest subfolder
    crest_dir_prefix = "crest_checks"
    dirs_lst = [dir for dir in os.listdir() if crest_dir_prefix in dir]
    folder_nums = []
    if not dirs_lst:
        crest_dir = crest_dir_prefix+"_1"
    else:
        for direc in dirs_lst:
            folder_nums.append(int(direc.split("_")[2]))
        crest_dir = f"{crest_dir_prefix}_{max(folder_nums) + 1}"
    os.system(f"mkdir -p {crest_dir}")
    print(f"\n####\nWorking in {crest_dir}\n####\n")

    crest_check = f'''cp {filename} {crest_dir}
                echo {spc_info[-2]} > {crest_dir}/.CHRG
                echo {int(spc_info[-1])-1} > {crest_dir}/.UHF
                cd {crest_dir}
                crest --for {filename} --prop singlepoint --ewin 50. &> crest_ouput.out
                '''
    with subprocess.Popen(crest_check, stdout=subprocess.PIPE, shell=True) as p:
        p.communicate()
        p.wait()

    with open(f"{crest_dir}/crest_ensemble.xyz","r",encoding="utf-8") as f:
        geo_list = from_xyz_trajectory_string(f.read())

    return [geo for geo,_ in geo_list]
