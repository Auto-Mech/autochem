""" some comparison functions
"""
import functools
import numpy
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates
from ._repr import coulomb_spectrum as _coulomb_spectrum


def almost_equal(geo1, geo2, rtol=2e-5):
    """ are these geometries numerically equal?
    """
    ret = False
    if _symbols(geo1) == _symbols(geo2):
        ret = numpy.allclose(_coordinates(geo1), _coordinates(geo2), rtol=rtol)
    return ret


def almost_equal_coulomb_spectrum(geo1, geo2, rtol=2e-5):
    """ do these geometries have similar coulomb spectrums?
    """
    ret = numpy.allclose(_coulomb_spectrum(geo1), _coulomb_spectrum(geo2),
                         rtol=rtol)
    return ret


def argunique_coulomb_spectrum(geos, seen_geos=(), rtol=2e-5):
    """ get indices of unique geometries, by coulomb spectrum
    """
    comp_ = functools.partial(almost_equal_coulomb_spectrum, rtol=rtol)
    idxs = _argunique(geos, comp_, seen_items=seen_geos)
    return idxs


def _argunique(items, comparison, seen_items=()):
    """ get the indices of unique items using some comparison function
    """
    idxs = []
    seen_items = list(seen_items)
    for idx, item in enumerate(items):
        if not any(comparison(item, seen_item) for seen_item in seen_items):
            idxs.append(idx)
            seen_items.append(item)
    idxs = tuple(idxs)
    return idxs
