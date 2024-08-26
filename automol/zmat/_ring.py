""" Level 4 Z-Matrix functions for generating ring information
"""

import math
import numpy
import os
import subprocess

from ..graph import base as graph_base
from ._conv import (
    distance,
    graph,
)
from .base import (
    coordinates,
    dihedral_angle_names,
    value_dictionary,
    from_geometry,
)

from ..geom.base import (
    dihedral_angle,
    from_xyz_trajectory_string,
    xyz_string,
    string,
)
from geom._ring import (
    translate_to_ring_center,
    get_displacement,
    ring_only_geometry,
)


# Get information for all rings at once
def all_rings_atoms(zma, tsg=None):
    """Get ring atoms.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    rings_atoms = graph_base.rings_atom_keys(graph(zma, dummy=True))
    if tsg is not None:
        for ring in graph_base.ts.forming_rings_bond_keys(tsg):
            # Determine number of atoms in the ring
            all_atoms = set()
            for ring_bnd in ring:
                all_atoms = all_atoms | ring_bnd
            natoms = len(all_atoms)

            # intialize list with indices from first bond
            atma, atmb = list(ring)[0]
            ring_atoms = [atma, atmb]

            # Iteratively add to ring idx list by finding with bnd has
            # the idx at end of current list to maintain connectivity
            while len(ring_atoms) != natoms:
                for ring_bnd in ring:
                    atma, atmb = ring_bnd
                    if atma == ring_atoms[-1] and atmb not in ring_atoms:
                        ring_atoms.append(atmb)
                    elif atmb == ring_atoms[-1] and atma not in ring_atoms:
                        ring_atoms.append(atma)

            # Check that ring is not already present in rings_atoms
            do_not_add_ring = 0
            for rng in rings_atoms:
                if not set(rng).difference(set(ring_atoms)):
                    do_not_add_ring = 1
            if do_not_add_ring:
                continue

            # Add to overall list
            rings_atoms = list(rings_atoms)
            # Added sort as I expect that connectivity is defined in "usual" way
            ring_atoms = sorted(ring_atoms)
            # Connectivity should still be preserved,
            # If it is not COME BACK HERE AND FIX

            # TODO - test
            # I think it is not ALWAYS preserved (see fused rings!!)
            #  so maybe I should remove the sort

            rings_atoms.append(tuple(ring_atoms))
            rings_atoms = frozenset(rings_atoms)

    return rings_atoms


def all_rings_distances(zma, rings_atoms):
    """For every ring present in the system. determine the
    distances between each pair of ring atoms.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """
    return tuple(ring_distances(zma, ring_atoms) for ring_atoms in rings_atoms)


def all_rings_distances_reasonable(zma, rings_atoms):
    """For every ring present in the system. determine the
    distances between each pair of ring atoms.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """
    # TODO HOW CAN IT EVER BE FALSE IF I CREATE VALUE DCT AND USE IT WITH SAME ZMA
    # Currently only used in tests
    condition = True
    for ring_atoms in rings_atoms:
        dist_val_dct = ring_distances(zma, ring_atoms)
        condition = ring_distances_reasonable(zma, ring_atoms, dist_val_dct)

    return condition


def all_rings_dihedrals(zma, rings_atoms):
    """Get ring dihedral names and their angle values

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """
    #return tuple(ring_dihedrals(zma, ring_atoms) for ring_atoms in rings_atoms)
    return tuple(complete_ring_dihedrals(zma, ring_atoms) for ring_atoms in rings_atoms)


def all_rings_dct(zma, rings_atoms):
    """Build a dictionary which relates the indices of the atoms
    of various rings to their dihedrals and sampling ranges.

    {rng_idx1-rng_idx2-rng_idx3: {Dn: [min, max], Dn2: [min, max]}}
    """

    ring_dct = {}
    for ring_atoms in rings_atoms:
        dct_label = "-".join(str(atm + 1) for atm in ring_atoms)
        ring_dct[dct_label] = ring_samp_ranges(zma, ring_atoms)

    return ring_dct


# Functions for a single ring
def ring_distances(zma, rng_atoms):
    """Return the distances between each pair of ring atoms.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    dist_value_dct = {}
    for i, _ in enumerate(rng_atoms):
        dist_value_dct[i] = distance(zma, rng_atoms[i - 1], rng_atoms[i])

    return dist_value_dct


def ring_distances_reasonable(zma, rng_atoms, dist_value_dct, thresh=0.3):
    """Are the distances between ring atoms reasonable?

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    condition = True
    for i, rng_atom in enumerate(rng_atoms):
        chk_dist = dist_value_dct[i] - distance(zma, rng_atoms[i - 1], rng_atom)
        if abs(chk_dist) > thresh:
            condition = False

    return condition


def ring_dihedrals(zma, rng_atoms):
    """Get ring dihedral names and their angle values

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    coos = coordinates(zma)
    da_names = dihedral_angle_names(zma)
    val_dct = value_dictionary(zma)

    ring_value_dct = {}
    for da_name in da_names:
        da_idxs = list(coos[da_name])[0]
        if len(list(set(da_idxs) & set(rng_atoms))) == 4:
            ring_value_dct[da_name] = val_dct[da_name]

    return ring_value_dct

def complete_ring_dihedrals(zma, rng_atoms):
    """Get ring dihedral names and their angle values
    of all the atoms involved in a ring, not only N-3 dihs

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    coos = coordinates(zma)
    da_names = dihedral_angle_names(zma)
    val_dct = value_dictionary(zma)

    ring_value_dct = {}
    for da_name in da_names:
        da_idxs = list(coos[da_name])[0]
        if da_idxs[0] in rng_atoms:
            ring_value_dct[da_name] = val_dct[da_name]

    return ring_value_dct


def ring_samp_ranges(zma, rng_atoms):
    """Set sampling range for ring dihedrals.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    samp_range_dct = {}
    ring_value_dct = ring_dihedrals(zma, rng_atoms)
    for key, value in ring_value_dct.items():
        samp_range_dct[key] = [value - math.pi / 4, value + math.pi / 4]

    return samp_range_dct


def checks_with_crest(filename,spc_info,vma,rings_atoms,eps=0.2):
    """Performs checks on ring geometries to determine unique sampling points for
    puckering protocol

    :param filename: input file containing geometries, molden style
    :type filename: string
    :param spc_info: input file containing geometries, molden style
    :type spc_info: automech data structure
    :param vma: value matrix from reference Z-Matrix
    :type vma: automol.vma object
    :param rings_atoms: list of all atoms in rings
    :type rings_atoms: list of lists of ints
    :output unique_zmas: Z-matrices of unique samples
    :type unique_zmas: list of automol.zmat objects
    """

    def z_score_normalize(features):
        """z-score normalization (0 mean, 1 std)

        :param features: array of features
        :type features: numpy.array
        :output normalized_features: array of normalized features
        :type normalized_features: numpy.array
        """
        mean = numpy.mean(features, axis=0)
        std = numpy.std(features, axis=0)
        normalized_features = (features - mean) / std
        return normalized_features

    def min_max_normalize(features):
        """min-max normalization (all values between 0 and 1)

        :param features: array of features
        :type features: numpy.array
        :output normalized_features: array of normalized features
        :type normalized_features: numpy.array
        """
        min_val = numpy.min(features, axis=0)
        max_val = numpy.max(features, axis=0)
        normalized_features = (features - min_val) / (max_val - min_val)
        return normalized_features

    def dbscan(features, eps, min_samples=1):
        """Density based clustering algorithm

        :param features: array of features
        :type features: numpy.array
        :param eps: pradius of a neighborhood centered on a given point
        :type eps: cluster
        :param min_samples: minimum number of points in cluster
        :type min_samples: int
        :output labels: list of labels for each data point
        :type labels: list of ints from 1 to n_clusters
        """

        def find_neighbors(features,point):
            distances = numpy.linalg.norm(features - features[point], axis=1)
            return numpy.asarray(distances <= eps).nonzero()[0]

        num_points = features.shape[0]
        labels = numpy.full(num_points, 0)  # Initialize all labels as 0 (unclassified)
        cluster_id = 1

        for point in range(num_points):
            #print(f"Working on point {point}, initial label is {labels[point]}")
            if labels[point] != 0:  # Skip if already classified
                continue

            # Find neighbors
            neighbors = find_neighbors(features,point)
            #print(f"Initial neighbors: {neighbors}")
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
                            neighbors = numpy.concatenate((neighbors, new_neighbors))
                    i += 1
                cluster_id += 1
            #print(f"Now label is {labels[point]}")
        print("labels",labels)
        return labels

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
                crest --for {filename} --prop singlepoint --ewin 100. &> crest_ouput.out
                '''
    with subprocess.Popen(crest_check, stdout=subprocess.PIPE, shell=True) as p:
        output, err = p.communicate()
        p_status = p.wait()

    with open(f"{crest_dir}/crest_ensemble.xyz","r",encoding="utf-8") as f:
        geo_list = from_xyz_trajectory_string(f.read())
    crest_geos = [geo for geo,_ in geo_list]
    print("rings_atoms: ", rings_atoms)
    sub_geos = [ring_only_geometry(geoi,rings_atoms) for geoi in crest_geos]

    subgeo_strings = []
    with open("sub_geoms.xyz","w",encoding="utf-8") as f:
        for geoi in sub_geos:
            geo_string = xyz_string(geoi, comment="  ")
            subgeo_strings.append(geo_string)
            f.write(geo_string+"\n")

    # Calculate RMSD for each pair of molecules
    import rdkit
    mols = [rdkit.Chem.rdmolfiles.MolFromXYZBlock(geoi) for geoi in subgeo_strings]
    num_mols = len(mols)
    rmsd_matrix = numpy.zeros((num_mols, num_mols))
    for i in range(num_mols):
        for j in range(i+1, num_mols):
            rmsd = rdkit.Chem.AllChem.GetBestRMS(mols[i], mols[j])
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd

    # Compute displacements from mean ring planes
    z_list = []
    for geoi in crest_geos:
        z_local = []
        for ring_atoms in rings_atoms:
            geo_string = string(geoi, angstrom=True)
            geo_list = [ [float(x) for x in line.split(
                        )[1:]] for line in geo_string.split('\n') ]
            coords = [xyz for i,xyz in enumerate(geo_list) if i in ring_atoms]
            coords = numpy.array(coords)
            coords = translate_to_ring_center(coords)
            z = get_displacement(coords)
            z_local.extend(z)
        z_list.append(z_local)
    input_z = numpy.array(z_list)

    # Compute dihedrals of rings atoms
    dih_list = []
    for geoi in crest_geos:
        dih_local = []
        for ring_atoms in rings_atoms:
            for i in range(len(ring_atoms)):
                idxs = [ring_atoms[j%len(ring_atoms)] for j in range(i,i+4)]
                dih = dihedral_angle(geoi,*idxs)
                if dih > numpy.pi:
                    dih -= 2*numpy.pi
                elif dih < -numpy.pi:
                    dih += 2*numpy.pi
                dih_local.append(dih)
        dih_list.append(dih_local)
    input_dih = numpy.array(dih_list)

    # Clustering with DBSCAN algorithm
    #input_z = min_max_normalize(input_z)
    features = numpy.hstack((input_dih,input_z,rmsd_matrix))
    labels = dbscan(features, eps=eps, min_samples=1)
    print("labels: ",labels)

    unique_geos = []
    visited_labels = set()
    for label,geoi in zip(labels,crest_geos):
        if label not in visited_labels:
            visited_labels.add(label)
            unique_geos.append(geoi)

    unique_zmas = [from_geometry(vma, geoi) for geoi in unique_geos]

    return unique_zmas
