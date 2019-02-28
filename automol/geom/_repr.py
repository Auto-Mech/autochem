""" molecular geometry representations
"""
import numpy
from more_itertools import unique_everseen as _unique_everseen
import phycon.elements as pce
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates


def formula(geo):
    """ molecular formula, as a dictionary
    """
    syms = list(map(pce.standard_case, _symbols(geo)))
    frm_dct = {sym: syms.count(sym) for sym in _unique_everseen(syms)}
    return frm_dct


def coulomb_spectrum(geo):
    """ (sorted) coulomb matrix eigenvalue spectrum
    """
    mat = _coulomb_matrix(geo)
    vals = tuple(sorted(numpy.linalg.eigvalsh(mat)))
    return vals


def _coulomb_matrix(geo):
    nums = numpy.array(list(map(pce.number, _symbols(geo))))
    xyzs = numpy.array(_coordinates(geo))

    _ = numpy.newaxis
    natms = len(nums)
    diag_idxs = numpy.diag_indices(natms)
    tril_idxs = numpy.tril_indices(natms, -1)
    triu_idxs = numpy.triu_indices(natms, 1)

    zxz = numpy.outer(nums, nums)
    rmr = numpy.linalg.norm(xyzs[:, _, :] - xyzs[_, :, :], axis=2)

    mat = numpy.zeros((natms, natms))
    mat[diag_idxs] = nums ** 2.4 / 2.
    mat[tril_idxs] = zxz[tril_idxs] / rmr[tril_idxs]
    mat[triu_idxs] = zxz[triu_idxs] / rmr[triu_idxs]
    return mat
