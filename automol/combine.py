"""
    Handle information for multiple species
      (1) van der Waals wells,
      (2) reaction path points
"""

import functools
import itertools
from operator import mul

from . import form, geom


# Combines via geometry
def formula_string(*geos):
    """get the overall combined stoichiometry"""
    fml = form.join_sequence(list(map(geom.formula, geos)))
    return form.string2(fml)


def fake_vdw_geometry(*geos):
    """put two geometries together in a fake well"""
    return geom.join_sequence(geos)


def fake_vdw_frequencies(*geos):
    """generate fake frequencies to fill in the missing internal modes of a fake well"""

    def _external_mode_count(geo):
        return 3 if geom.is_atom(geo) else 5 if geom.is_linear(geo) else 6

    fake_geo = fake_vdw_geometry(*geos)
    nfake = sum(map(_external_mode_count, geos)) - _external_mode_count(fake_geo)
    return tuple(10.0 * (i + 1) for i in range(nfake))


# Combines by elec levels
def electronic_energy_levels(*elec_levels_lst):
    """Combine electronic levels for multiple species"""

    def _merge_level_combinations(elec_level_comb):
        """Merge a combination of fragment electronic levels into one"""
        enes, mults = zip(*elec_level_comb)
        return (sum(enes), functools.reduce(mul, mults))

    # 1. Get combinations
    elec_level_combs = list(itertools.product(*elec_levels_lst))
    # 2. Combine them
    elec_levels = tuple(map(_merge_level_combinations, elec_level_combs))
    # 3. Sort them
    elec_levels = tuple(sorted(elec_levels, key=lambda x: x[0]))
    # 4. Combine duplicates
    # Note: Must account for the increased degeneracy!
    elec_levels = tuple(
        (e, m * len(list(g))) for (e, m), g in itertools.groupby(elec_levels)
    )

    return elec_levels
