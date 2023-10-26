"""A data structure for describing a multi-dimensional potential

Values along the potential are stored in an xarray.Dataset object. Energies are required
for every potential object, but it can additionally store geometry, gradient, hessian,
or z-matrix values along the potential.

This allows accessing values along the potential using the coordinate values.

For example,

    >>> ene = pot["energy"].sel(q1=0.1, q2=0.25, method="nearest")
    >>> ene = potent.value(pot, 0.1, 0.25)   # using module function

would grab the energy nearest to the coordinates (q1, q2) = (0.1, 0.25), whereas

    >>> geo = pot["geom"].sel(q1=0.1, q2=0.25, method="nearest")
    >>> ene = potent.value(pot, 0.1, 0.25, key="geom")   # using module function

would grab the geometry for the same point.

To scale the energies in the potential by 5, one could simply do

    >>> pot["energy"] *= 5
    >>> pot = potent.scaled(pot, 5)   # using module function

which would scale the energies in place.
"""
from typing import Dict, List, Optional

import numpy
import xarray

Potential = xarray.Dataset

ALLOWED_KEYS = {"geom", "grad", "hess", "zmat"}
KEY_TYPE_DCT = {
    "energy": float,
    # eventually, these could be replaced with specific types
    "geom": object,
    "grad": object,
    "hess": object,
    "zmat": object,
}


# Constructors
def from_dict(
    ene_dct: Dict[tuple, float],
    coord_names: Optional[List[str]] = None,
    aux_dct_dct: Optional[Dict[str, Dict[tuple, object]]] = None,
) -> Potential:
    """Construct a potential data structure from a dictionary of energy values

    :param ene_dct: A dictionary of energy values, by coordinate value
    :type ene_dct: Dict[tuple, float]
    :param coord_names: Optionally, specify the names for each coordinate in order
    :type coord_names: List[str], optional
    :param aux_dct_dct: Dictionaries of auxiliary values along the potential:
        {
            "geom": {
                (0.15, 0.1): <automol geom data structure>,
                (0.30, 0.2): <automol geom data structure>,
                ...
            }, ...
        }
        the keywords/types accepted are: "geom", "grad", "hess", and "zmat"
    :type aux_dct_dct: Optional[Dict[str, Dict[tuple, object]]]
    :returns: A potential data structure
    :rtype: Potential
    """
    coord_vals_lst = list(map(sorted, map(set, zip(*ene_dct.keys()))))
    shape_ = tuple(map(len, coord_vals_lst))

    def array_from_dict_(val_dct, dtype):
        """Create an array from a dictionary of values along the potential"""
        val_arr = numpy.empty(shape_, dtype=dtype)
        for idxs, _ in numpy.ndenumerate(val_arr):
            coords = tuple(c[i] for i, c in zip(idxs, coord_vals_lst))
            val_arr[idxs] = val_dct[coords]
        return val_arr

    ene_arr = array_from_dict_(ene_dct, float)

    if aux_dct_dct is None:
        aux_arr_dct = None
    else:
        aux_arr_dct = {k: array_from_dict_(v, object) for k, v in aux_dct_dct.items()}

    return from_data(
        ene_arr=ene_arr,
        coord_vals_lst=coord_vals_lst,
        coord_names=coord_names,
        aux_arr_dct=aux_arr_dct,
    )


def from_data(
    ene_arr: List[List[float]],
    coord_vals_lst: List[List[float]],
    coord_names: Optional[List[str]] = None,
    aux_arr_dct: Optional[List[List[object]]] = None,
) -> Potential:
    """Construct a potential data structure from coordinate and energy values

    :param ene_arr: The array of energy values for the potential, with shape matching
        the numbers of coordinates in order
    :type pot_vals: List[List[float]]
    :param coord_vals_lst: Lists of values for each coordinate
    :type coord_vals_lst: List[List[float]]
    :param coord_names: Optionally, specify the names for each coordinate
    :type coord_names: List[str], optional
    :param aux_arr_dct: A dictionary of arrays of auxiliary values along the potential;
        the allowed keys are: "geom", "grad", "hess", and "zmat"
    :returns: A potential data structure
    :rtype: Potential
    """
    shape_ = tuple(map(len, coord_vals_lst))

    coord_names = (
        [f"q{i}" for i, _ in enumerate(coord_vals_lst)]
        if coord_names is None
        else coord_names
    )

    def array_(vals, dtype):
        arr = numpy.empty(shape_, dtype=dtype)
        arr[:] = vals
        return arr

    data_vars = {"energy": (coord_names, array_(ene_arr, float))}

    if aux_arr_dct is not None:
        keys_ = set(aux_arr_dct)
        assert (
            keys_ <= ALLOWED_KEYS
        ), f"Non-allowed keys detected in aux values: {keys_} !<= {ALLOWED_KEYS}"

        data_vars.update(
            {k: (coord_names, array_(v, object)) for k, v in aux_arr_dct.items()}
        )

    return xarray.Dataset(
        data_vars=data_vars, coords=dict(zip(coord_names, coord_vals_lst))
    )


# Conversions
def dict_(
    pot: Potential, key: str = "energy", index: bool = False, drop_null: bool = False
) -> Dict[tuple, float]:
    """Convert a potential to a dictionary

    :param pot: A potential
    :type pot: Potential
    :param key: Which values to select, defaults to "energy"
    :type key: str, optional
    :param index: Use indices instead of coordinates for keys?, defaults to False
    :type index: bool, optional
    :param drop_null: Drop null values?
    :type drop_null: bool, optional
    :return: The potential formatted as a dictionary
    :rtype: Dict[tuple, float]
    """

    def is_null_(val):
        if KEY_TYPE_DCT[key] is float:
            return numpy.isnan(val)
        return val is None

    if key not in keys(pot):
        return None

    coord_vals_lst = coordinates_values(pot)

    pot_dct = {}
    pot_vals = values(pot, key=key)
    for idxs, val in numpy.ndenumerate(pot_vals):
        key_ = idxs if index else tuple(c[i] for i, c in zip(idxs, coord_vals_lst))
        if not is_null_(val) or not drop_null:
            pot_dct[key_] = val

    return pot_dct


# Getters
def shape(pot: Potential) -> List[int]:
    """Get the shape of a Potential object

    :param pot: A potential
    :type pot: Potential
    :return: The shape
    :rtype: List[int]
    """
    return pot.energy.shape


def keys(pot: Potential) -> List[str]:
    """Get a list of keys for the values stored in a Potential object

    :param pot: A potential
    :type pot: Potential
    :return: The potential keys
    :rtype: List[str]
    """
    return tuple(pot.keys())


def values(pot: Potential, key: str = "energy") -> numpy.ndarray:
    """Get the potential values in a Potential object

    :param pot: A potential
    :type pot: Potential
    :param key: Which values to select, defaults to "energy"
    :type key: str, optional
    :return: The potential values
    :rtype: numpy.ndarray
    """
    if key not in keys(pot):
        return None

    return pot[key].values


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
    return pot.energy.dims


def value(
    pot: Potential,
    *coords: float,
    key: str = "energy",
    method: str = "nearest",
    **coord_dct: float,
) -> float:
    """Get the value of the potential at a specific set of coordinates

    Coordinates can be passed in either directly or using their coordinate names:

        value(pot, 0.25, 0.333)
        value(pot, q1=0.25, q2=0.333)

    :param pot: A potential
    :type pot: Potential
    :param coords: The coordinates to get the potential value for
    :type coords: List[float]
    :param key: Which values to select, defaults to "energy"
    :type key: str, optional
    :param method: The method of value selection, defaults to "nearest"
        (See xarray.DataArray.sel for other options)
    :type method: str, optional
    :param coord_dct: The coordinates to get the potential value for
    :type coord_dct: Dict[str, float]
    :return: The potential value for these coordinates
    :rtype: float
    """
    if key not in keys(pot):
        return None

    assert not (
        coords and coord_dct
    ), f"Coordinates must be given by value *or* by keyword:\n{coords}\n{coord_dct}"

    coord_dct = dict(zip(coordinate_names(pot), coords)) if not coord_dct else coord_dct
    return pot[key].sel(**coord_dct, method=method).item()


# Value checking and equality
def almost_equal(
    pot1: Potential, pot2: Potential, rtol: float = 2e-3, atol: float = 2e-6
) -> bool:
    """Are these two potentials almost equal?

    Currently not comparing auxiliary values

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
    pot_vals1, pot_vals2 = map(values, [pot1, pot2])
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
    return bool(pot.energy.notnull().any())


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
    pot = pot.copy()  # could do deep=True, but only the energies need to be copied
    pot["energy"] = pot.energy.copy()
    pot["energy"] *= factor
    return pot


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
    pot = pot.copy()  # could do deep=True, but only the energies need to be copied
    pot["energy"] = pot.energy.copy()
    pot["energy"] *= 1 / (1 + squash_factor * pot["energy"])
    return pot
