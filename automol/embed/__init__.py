""" geometry embedding using the distance geometry algorithm

This module does the numerical work behind the functions in
automol.graph.embed, which allow one to determine a qualitatively appropriate
geometry for a molecular graph.

The algorithm used is the "Distance Geomtry" Algorithm, which generates
approximate geometries by a randomized guess at the distance matrix, within
heuristic distance bounds, which is then converted to an approximate geometry
that satifies these distances.

See Blaney, J. M.; Dixon, J. S. “Distance Geometry in Molecular Modeling”.
Reviews in Computational Chemistry; VCH: New York, 1994.

The steps in the algorithm are as follows:

    1. Generate distance bounds matrix B. (Here, B is replaced with L and U).
    2. Narrow the bounds in B by triangle smoothing.
    3. Generate a distance matrix D by uniform sampling within the bounds.
    4. Generate the metric matrix G (matrix of position vector dot products).
    5. Diagonalize G and determine the principal components.
    6. The three largest eigenvectors and eigenvalues of G can be used to
    generate x, y, z coordinates for the molecule which approximately
    correspond to the distance matrix D.
    7. Do error refinement to clean up the structure and enforce correct
    chirality.

Step 1. is done directly from the graph in automol.graph.embed.

Steps 2-6. constitute the core of the distance geometry algorithm and are
implemented here in _dgeom.py.

Step 7. cleans up the result of the distance geometry and is done in
_cleanup.py.
"""
from automol.embed._dgeom import triangle_smooth_bounds_matrices
from automol.embed._dgeom import sample_distance_matrix
from automol.embed._dgeom import distances_from_center
from automol.embed._dgeom import metric_matrix
from automol.embed._dgeom import coordinates_from_metric_matrix
from automol.embed._dgeom import metric_matrix_from_coordinates
from automol.embed._dgeom import distance_matrix_from_coordinates


__all__ = [
    'triangle_smooth_bounds_matrices',
    'sample_distance_matrix',
    'distances_from_center',
    'metric_matrix',
    'coordinates_from_metric_matrix',
    'metric_matrix_from_coordinates',
    'distance_matrix_from_coordinates',
]
