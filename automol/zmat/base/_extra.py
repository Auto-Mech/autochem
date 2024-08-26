""" Level 2 functions for sampling and related actions
"""

import os
import subprocess
import itertools

import numpy

from . import coordinates as zcoords

from ._core import (
    set_values_by_name,
    value_dictionary,
    from_geometry,
)
from ._core import value, set_key_matrix
from ...geom import dihedral_angle
#from ...zmat import coordinates as zcoords
from ...vmat import coordinates, key_matrix

from ...geom.base import (
    from_xyz_trajectory_string,
    xyz_string,
    string,
)

# from ...geom._ring import (
#     translate_to_ring_center,
#     get_displacement,
#     ring_only_geometry,
# )

def samples(zma, nsamp, range_dct):
    """randomly sample over torsional dihedrals"""
    _names = tuple(range_dct.keys())
    ranges = tuple(range_dct.values())
    vals_lst = _sample_over_ranges(ranges, nsamp)

    zmas = tuple(
        set_values_by_name(zma, dict(zip(_names, vals)), degree=False)
        for vals in vals_lst
    )

    return zmas


def samples_avg_dih(zma, geo, tors_dcts, average_dih,ring_tors_dct,dih_remover):
    """randomly sample over torsional dihedrals"""
    vals_lst = []
    iter_combos = []
    _names, avg_dih = [], []

    all_sampled_torsions = []
    for _,samp_range_dct in tors_dcts:
        all_sampled_torsions.extend(list(samp_range_dct.keys()))

    for key_dct,samp_range_dct in tors_dcts:

        repeats = len(samp_range_dct.keys())
        ring_atoms = [int(idx)-1 for idx in key_dct.split('-')]

        if repeats != len(ring_atoms):
            keymat = [list(item) for item in(key_matrix(zma))]
            # Find all DHs of ring atoms from 4th atom
            # is DH in tors-dcts?
            changed_dh = []

            for at in ring_atoms[3:]:
                dih=f"D{at}"
                if dih not in [b for a,b in dih_remover if a != key_dct]:
                    if dih not in all_sampled_torsions:
                        changed_dh.append(dih)
                        idx = ring_atoms.index(at)
                        keymat[at] = [ring_atoms[idx-1], \
                                      ring_atoms[idx-2], ring_atoms[idx-3]]

                        ring_tors_dct[key_dct].update({dih:(0.,0.)})

            zma = set_key_matrix(zma, keymat)

            # rebuild the zmatrix with new value of that dih
            key_coord_dct = zcoords(zma)
            new_key_dct = {}
            for name,coos in key_coord_dct.items():
                atm_idxs = coos[0]
                if len(atm_idxs) == 2:
                    new_key_dct[name] = value(zma, name, angstrom=True)
                elif len(atm_idxs) == 3:
                    new_key_dct[name] = value(zma, name, degree=True)
                elif len(atm_idxs) == 4:
                    new_key_dct[name] = value(zma, name, degree=True)
                    if name in changed_dh: # compute DH with previous three ring atoms
                        new_key_dct[name] = dihedral_angle(geo, *atm_idxs,degree=True)
            zma = set_values_by_name(zma, new_key_dct)


    for key_dct,samp_range_dct in tors_dcts:

        repeats = len(samp_range_dct.keys())

        _names.extend(list(samp_range_dct.keys()))
        iter_combos.append(itertools.product([1,0,-1],repeat=repeats))
        avg_dih.append(average_dih[key_dct])

    samp_mat = itertools.product(*iter_combos)

    for samp in samp_mat:
        vals = []
        for i,dih_value in enumerate(avg_dih):
            vals.extend( [val * dih_value for val in samp[i]] )
        vals_lst.append(tuple(vals))

    zmas = tuple(
        set_values_by_name(zma, dict(zip(_names, vals)), degree=False)
        for vals in vals_lst
    )

    return zmas


def constraint_dict(zma, const_names, var_names=()):
    """Build a dictionary of constraints

    has to round the values for the filesystem
    """

    # Get the list names sorted for dictionary
    rnames = (name for name in const_names if "R" in name)
    anames = (name for name in const_names if "A" in name)
    dnames = (name for name in const_names if "D" in name)
    rnames = tuple(sorted(rnames, key=lambda x: int(x.split("R")[1])))
    anames = tuple(sorted(anames, key=lambda x: int(x.split("A")[1])))
    dnames = tuple(sorted(dnames, key=lambda x: int(x.split("D")[1])))
    constraint_names = rnames + anames + dnames

    # Remove the scan coordinates so they are not placed in the dict
    constraint_names = tuple(name for name in constraint_names if name not in var_names)

    # Build dictionary
    if constraint_names:
        zma_vals = value_dictionary(zma)
        zma_coords = coordinates(zma)
        assert set(constraint_names) <= set(zma_coords.keys()), (
            "Attempting to constrain coordinates not in zma:"
            f"\n{constraint_names}\n{zma_coords}"
        )
        _dct = dict(
            zip(
                constraint_names,
                (round(zma_vals[name], 2) for name in constraint_names),
            )
        )
    else:
        _dct = None

    return _dct


def set_constraint_names(zma, tors_names, tors_model):
    """Determine the names of constraints along a torsion scan"""

    constraint_models = ("1dhrf", "1dhrfa", "tau-1dhrf", "tau-1dhrfa")

    const_names = tuple()
    if tors_names and tors_model in constraint_models:
        if tors_model in ("1dhrf", "tau-1dhrf"):
            const_names = tuple(itertools.chain(*tors_names))
        elif tors_model in ("1dhrfa", "tau-1dhrfa"):
            coords = list(coordinates(zma))
            const_names = tuple(coord for coord in coords)

    return const_names


# helpers
def _sample_over_ranges(rngs, nsamp):
    """randomly sample over several ranges"""
    nrng = len(rngs)
    samp_mat = numpy.random.rand(nsamp, nrng)
    for i, (start, stop) in enumerate(rngs):
        samp_mat[:, i] = samp_mat[:, i] * (stop - start) + start
    return tuple(map(tuple, samp_mat))


# def checks_with_crest(filename,spc_info,vma,rings_atoms,eps=0.2):
#     """Performs checks on ring geometries to determine unique sampling points for
#     puckering protocol

#     :param filename: input file containing geometries, molden style
#     :type filename: string
#     :param spc_info: input file containing geometries, molden style
#     :type spc_info: automech data structure
#     :param vma: value matrix from reference Z-Matrix
#     :type vma: automol.vma object
#     :param rings_atoms: list of all atoms in rings
#     :type rings_atoms: list of lists of ints
#     :output unique_zmas: Z-matrices of unique samples
#     :type unique_zmas: list of automol.zmat objects
#     """

#     def z_score_normalize(features):
#         """z-score normalization (0 mean, 1 std)

#         :param features: array of features
#         :type features: numpy.array
#         :output normalized_features: array of normalized features
#         :type normalized_features: numpy.array
#         """
#         mean = numpy.mean(features, axis=0)
#         std = numpy.std(features, axis=0)
#         normalized_features = (features - mean) / std
#         return normalized_features

#     def min_max_normalize(features):
#         """min-max normalization (all values between 0 and 1)

#         :param features: array of features
#         :type features: numpy.array
#         :output normalized_features: array of normalized features
#         :type normalized_features: numpy.array
#         """
#         min_val = numpy.min(features, axis=0)
#         max_val = numpy.max(features, axis=0)
#         normalized_features = (features - min_val) / (max_val - min_val)
#         return normalized_features

#     def dbscan(features, eps, min_samples=1):
#         """Density based clustering algorithm

#         :param features: array of features
#         :type features: numpy.array
#         :param eps: pradius of a neighborhood centered on a given point
#         :type eps: cluster
#         :param min_samples: minimum number of points in cluster
#         :type min_samples: int
#         :output labels: list of labels for each data point
#         :type labels: list of ints from 1 to n_clusters
#         """

#         def find_neighbors(features,point):
#             distances = numpy.linalg.norm(features - features[point], axis=1)
#             return numpy.asarray(distances <= eps).nonzero()[0]

#         num_points = features.shape[0]
#         labels = numpy.full(num_points, 0)  # Initialize all labels as 0 (unclassified)
#         cluster_id = 1

#         for point in range(num_points):
#             #print(f"Working on point {point}, initial label is {labels[point]}")
#             if labels[point] != 0:  # Skip if already classified
#                 continue

#             # Find neighbors
#             neighbors = find_neighbors(features,point)
#             #print(f"Initial neighbors: {neighbors}")
#             if len(neighbors) < min_samples:
#                 labels[point] = -1  # Mark as noise
#             else:
#                 labels[point] = cluster_id
#                 i = 0
#                 while i < len(neighbors):
#                     neighbor_point = neighbors[i]
#                     if labels[neighbor_point] == -1:  # Noise point in cluster
#                         labels[neighbor_point] = cluster_id
#                     elif labels[neighbor_point] == 0:  # New unvisited point
#                         labels[neighbor_point] = cluster_id
#                         # Find more neighbors of this point
#                         new_neighbors = find_neighbors(features,neighbor_point)
#                         if len(new_neighbors) >= min_samples:
#                             neighbors = numpy.concatenate((neighbors, new_neighbors))
#                     i += 1
#                 cluster_id += 1
#             #print(f"Now label is {labels[point]}")
#         print("labels",labels)
#         return labels

#     # Setup crest subfolder
#     crest_dir_prefix = "crest_checks"
#     dirs_lst = [dir for dir in os.listdir() if crest_dir_prefix in dir]
#     folder_nums = []
#     if not dirs_lst:
#         crest_dir = crest_dir_prefix+"_1"
#     else:
#         for direc in dirs_lst:
#             folder_nums.append(int(direc.split("_")[2]))
#         crest_dir = f"{crest_dir_prefix}_{max(folder_nums) + 1}"
#     os.system(f"mkdir -p {crest_dir}")
#     print(f"\n####\nWorking in {crest_dir}\n####\n")

#     crest_check = f'''cp {filename} {crest_dir}
#                 echo {spc_info[-2]} > {crest_dir}/.CHRG
#                 echo {int(spc_info[-1])-1} > {crest_dir}/.UHF
#                 cd {crest_dir}
#                 crest --for {filename} --prop singlepoint --ewin 100. &> crest_ouput.out
#                 '''
#     with subprocess.Popen(crest_check, stdout=subprocess.PIPE, shell=True) as p:
#         output, err = p.communicate()
#         p_status = p.wait()

#     with open(f"{crest_dir}/crest_ensemble.xyz","r",encoding="utf-8") as f:
#         geo_list = from_xyz_trajectory_string(f.read())
#     crest_geos = [geo for geo,_ in geo_list]
#     print("rings_atoms: ", rings_atoms)
#     sub_geos = [ring_only_geometry(geoi,rings_atoms) for geoi in crest_geos]

#     subgeo_strings = []
#     with open("sub_geoms.xyz","w",encoding="utf-8") as f:
#         for geoi in sub_geos:
#             geo_string = xyz_string(geoi, comment="  ")
#             subgeo_strings.append(geo_string)
#             f.write(geo_string+"\n")

#     # Calculate RMSD for each pair of molecules
#     import rdkit
#     mols = [rdkit.Chem.rdmolfiles.MolFromXYZBlock(geoi) for geoi in subgeo_strings]
#     num_mols = len(mols)
#     rmsd_matrix = numpy.zeros((num_mols, num_mols))
#     for i in range(num_mols):
#         for j in range(i+1, num_mols):
#             rmsd = rdkit.Chem.AllChem.GetBestRMS(mols[i], mols[j])
#             rmsd_matrix[i, j] = rmsd
#             rmsd_matrix[j, i] = rmsd

#     # Compute displacements from mean ring planes
#     z_list = []
#     for geoi in crest_geos:
#         z_local = []
#         for ring_atoms in rings_atoms:
#             geo_string = string(geoi, angstrom=True)
#             geo_list = [ [float(x) for x in line.split(
#                         )[1:]] for line in geo_string.split('\n') ]
#             coords = [xyz for i,xyz in enumerate(geo_list) if i in ring_atoms]
#             coords = numpy.array(coords)
#             coords = translate_to_ring_center(coords)
#             z = get_displacement(coords)
#             z_local.extend(z)
#         z_list.append(z_local)
#     input_z = numpy.array(z_list)

#     # Compute dihedrals of rings atoms
#     dih_list = []
#     for geoi in crest_geos:
#         dih_local = []
#         for ring_atoms in rings_atoms:
#             for i in range(len(ring_atoms)):
#                 idxs = [ring_atoms[j%len(ring_atoms)] for j in range(i,i+4)]
#                 dih = dihedral_angle(geoi,*idxs)
#                 if dih > numpy.pi:
#                     dih -= 2*numpy.pi
#                 elif dih < -numpy.pi:
#                     dih += 2*numpy.pi
#                 dih_local.append(dih)
#         dih_list.append(dih_local)
#     input_dih = numpy.array(dih_list)

#     # Clustering with DBSCAN algorithm
#     #input_z = min_max_normalize(input_z)
#     features = numpy.hstack((input_dih,input_z,rmsd_matrix))
#     labels = dbscan(features, eps=eps, min_samples=1)
#     print("labels: ",labels)

#     unique_geos = []
#     visited_labels = set()
#     for label,geoi in zip(labels,crest_geos):
#         if label not in visited_labels:
#             visited_labels.add(label)
#             unique_geos.append(geoi)

#     unique_zmas = [from_geometry(vma, geoi) for geoi in unique_geos]

#     return unique_zmas
