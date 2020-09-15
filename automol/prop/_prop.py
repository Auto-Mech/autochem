"""
  Functions to handle
"""

import numpy


def total_dipole_moment(vec):
    """ Calculate the total dipole_moment value from the vector. Value
        is simply the norm of the vector.

        :param vector: dipole_moment vector in XYZ coords (Bohr)
        :type vector: tuple
        :rtype: float
    """

    assert len(vec) == 3

    return numpy.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)


def total_polarizability(tensor):
    """ Calculate the total static polarizability value from the tensor
        which is given in x,y,z coordinates (Bohr).

        Assumes that value is an average of the xx, yy, zz components.

        :param tensor: 3x3 polarizability tensor in XYZ coords (Bohr)
        :type tensor: tuple
        :rtype: float
    """

    assert len(tensor) == 3
    assert all(len(row) == 3 for row in tensor)

    return (1.0/3.0) * (tensor[0][0] + tensor[1][1] + tensor[2][2])
