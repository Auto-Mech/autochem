"""A data structure for describing a multi-dimensional potential
"""
from typing import Dict, List, Optional

import numpy
import xarray

Potential = xarray.Dataset


# Constructors
def from_dict(
    pot_dct: Dict[tuple, float],
    coord_names: Optional[List[str]] = None,
    geo_dct: Optional[Dict[tuple, object]] = None,
) -> Potential:
    """Construct a potential data structure from a dictionary of potential values

    :param pot_dct: A dictionary of potential values, by coordinate value
    :type pot_dct: Dict[tuple, float]
    :param pot_vals: The potential values, as an array with shape matching the numbers
        of coordinates in order
    :type pot_vals: List
    :param geo_dct: A dictionary of geometries along the potential, by coordinate vlue
    :type geo_dct: Dict[tuple, automol geom data structure]
    :returns: A potential data structure
    :rtype: Potential
    """
    coord_vals_lst = list(map(sorted, map(set, zip(*pot_dct.keys()))))

    shape_ = tuple(map(len, coord_vals_lst))
    pot_arr = numpy.array(list(pot_dct.values()), dtype=float)
    pot_arr = numpy.reshape(pot_arr, shape_)

    if geo_dct is None:
        geo_arr = None
    else:
        geo_arr = numpy.empty_like(pot_arr, dtype=object)
        idxs_lst, _ = zip(*numpy.ndenumerate(pot_arr))
        for idxs, geo in zip(idxs_lst, geo_dct.values()):
            geo_arr[idxs] = geo

    return from_data(
        pot_arr=pot_arr,
        coord_vals_lst=coord_vals_lst,
        coord_names=coord_names,
        geo_arr=geo_arr,
    )


def from_data(
    pot_arr: List[List[float]],
    coord_vals_lst: List[List[float]],
    coord_names: Optional[List[str]] = None,
    geo_arr: Optional[List[List[object]]] = None,
) -> Potential:
    """Construct a potential data structure from coordinate and potential values

    :param pot_vals: The potential values, as an array with shape matching the numbers
        of coordinates in order, or the corresponding flattened array
    :type pot_vals: List[List[float]]
    :param coord_vals_lst: Lists of values for each coordinate
    :type coord_vals_lst: List[List[float]]
    :param coord_names: Optionally, specify the names for each coordinate
    :type coord_names: List[str], optional
    :returns: A potential data structure
    :rtype: Potential
    """
    shape_ = tuple(map(len, coord_vals_lst))
    pot_arr = numpy.array(pot_arr, dtype=float)

    if coord_names is None:
        coord_names = [f"q{i}" for i, _ in enumerate(coord_vals_lst)]

    assert (
        numpy.shape(pot_arr) == shape_
    ), f"Potential values shape don't match coordinates:\n{pot_arr}\n{coord_vals_lst}"
    assert len(coord_names) == len(coord_vals_lst)

    data_vars = {"potential": (coord_names, pot_arr)}

    if geo_arr is not None:
        geo_arr_data = geo_arr
        geo_arr = numpy.empty_like(pot_arr, dtype=object)
        geo_arr[:] = geo_arr_data

        data_vars["geometry"] = (coord_names, geo_arr)

    return xarray.Dataset(
        data_vars=data_vars, coords=dict(zip(coord_names, coord_vals_lst))
    )


# Conversions
def dict_(
    pot: Potential, index: bool = False, drop_null: bool = False
) -> Dict[tuple, float]:
    """Convert a potential to a dictionary

    :param pot: A potential
    :type pot: Potential
    :param index: Use indices instead of coordinates for keys?, defaults to False
    :type index: bool, optional
    :param drop_null: Drop null values?
    :type drop_null: bool, optional
    :return: The potential formatted as a dictionary
    :rtype: Dict[tuple, float]
    """
    coord_vals_lst = coordinates_values(pot)

    pot_dct = {}
    pot_vals = potential_values(pot)
    for idxs, val in numpy.ndenumerate(pot_vals):
        key = idxs if index else tuple(c[i] for i, c in zip(idxs, coord_vals_lst))
        if not numpy.isnan(val) or not drop_null:
            pot_dct[key] = val

    return pot_dct


# Getters
def shape(pot: Potential) -> List[int]:
    """Get the shape of a Potential object

    :param pot: A potential
    :type pot: Potential
    :return: The shape
    :rtype: List[int]
    """
    return pot.potential.shape


def potential_values(pot: Potential) -> numpy.ndarray:
    """Get the potential values in a Potential object

    :param pot: A potential
    :type pot: Potential
    :return: The potential values
    :rtype: numpy.ndarray
    """
    return pot.potential.values


def coordinates_values(pot: Potential) -> List[numpy.ndarray]:
    """Get the values for each coordinate in a Potential object

    :param pot: A potential
    :type pot: Potential
    :return: The values for each coordinate
    :rtype: List[numpy.array]
    """
    return tuple(v.to_numpy() for _, v in pot.coords.items())


def coordinate_names(pot: Potential) -> List[str]:
    """Get the name of each coordinate in a Potential object

    :param pot: A potential
    :type pot: Potential
    :return: The name of each coordinate
    :rtype: List[numpy.array]
    """
    return pot.potential.dims


def value(
    pot: Potential, *coords: float, method: str = "nearest", **coord_dct: float
) -> float:
    """Get the value of the potential at a specific set of coordinates

    Coordinates can be passed in either directly or using their coordinate names:

        value(pot, 0.25, 0.333)
        value(pot, q1=0.25, q2=0.333)

    :param pot: A potential
    :type pot: Potential
    :param coords: The coordinates to get the potential value for
    :type coords: List[float]
    :param method: The method of value selection, defaults to "nearest"
        (See xarray.DataArray.sel for other options)
    :type method: str, optional
    :param coord_dct: The coordinates to get the potential value for
    :type coord_dct: Dict[str, float]
    :return: The potential value for these coordinates
    :rtype: float
    """
    assert not (
        coords and coord_dct
    ), f"Coordinates must be given by value *or* by keyword:\n{coords}\n{coord_dct}"

    coord_dct = dict(zip(coordinate_names(pot), coords)) if not coord_dct else coord_dct
    return float(pot.sel(**coord_dct, method=method).potential)


# Value checking and equality
def almost_equal(pot1: Potential, pot2: Potential, rtol=2e-3, atol=2e-6) -> bool:
    """Are these two potentials almost equal?

    :param pot1: A potential
    :type pot1: Potential
    :param pot2: Another potential
    :type pot2: Potential
    :param rtol: Relative tolerance for testing equality, defaults to 2e-3
    :type rtol: float, optional
    :param atol: Absolute tolerance for testing equality, defaults to 2e-6
    :type atol: float, optional
    :return: `True` if they are, `False` if they aren't
    :rtype: bool
    """
    # Check the shape
    shape1, shape2 = map(shape, [pot1, pot2])
    if shape1 != shape2:
        return False

    # Check the coordinate names
    coord_names1, coord_names2 = map(coordinate_names, [pot1, pot2])
    if coord_names1 != coord_names2:
        return False

    # Check the potential values
    pot_vals1, pot_vals2 = map(potential_values, [pot1, pot2])
    if not numpy.allclose(pot_vals1, pot_vals2, rtol=rtol, atol=atol):
        return False

    # Check the coordinate values
    coord_vals_lst1, coord_vals_lst2 = map(coordinates_values, [pot1, pot2])
    for coord_vals1, coord_vals2 in zip(coord_vals_lst1, coord_vals_lst2):
        if not numpy.allclose(coord_vals1, coord_vals2, rtol=rtol, atol=atol):
            return False

    return True


def has_defined_values(pot: Potential) -> bool:
    """Does this potential have non-null values?

    :param pot: A potential
    :type pot: Potential
    :return: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    return bool(pot.potential.notnull().any())


# Transformations
def scaled(pot: Potential, factor: float) -> Potential:
    """Scale a potential by some factor

    :param pot: A potential
    :type pot: Potential
    :param factor: The scale factor
    :type factor: float
    :return: The scaled potential
    :rtype: Potential
    """
    return pot * factor


def squashed(pot: Potential, squash_factor: float = 0.07) -> Potential:
    """Squash the potential according to the formula V <- V / (1 + f * V), where `f` is
    the "squash factor"

    :param pot: A potential
    :type pot: Potential
    :param factor: The scale factor, defaults to 0.07
    :type factor: float, optional
    :return: The scaled potential
    :rtype: Potential
    """
    return pot / (1 + squash_factor * pot)
