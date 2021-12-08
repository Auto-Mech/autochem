""" canonicalization functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!

Reference:
Schneider, Sayle, Landrum. J. Chem. Inf. Model. 2015, 55, 10, 2111–2120
"""
import itertools
import numpy
from automol.util import dict_
from automol.graph.base._core import atom_keys
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import implicit
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import atomic_numbers
from automol.graph.base._core import mass_numbers
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atoms_bond_keys
from automol.graph.base._algo import is_connected


# def canonical_keys(gra, backbone_only=True):
def canonical_keys(gra):
    """ Determine canonical keys for this graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param backbone_only: Consider backbone atoms only?
        :type backbone_only: bool
        :rtype: dict[int: int]
    """
    assert is_connected(gra), "Cannot canonicalize disconnected graph."

    gra = without_dummy_atoms(gra)
    gra = implicit(gra)

    # 1. Initial partition: all atoms in one class, with index 0.
    # idx_dct maps individual atom keys onto their class indices
    idx_dct = dict_.by_key({}, atom_keys(gra), fill_val=0)

    # 2. Relax class indices based on atom invariants
    idx_dct = relax_class_indices(
        gra, idx_dct, srt_eval_=sort_evaluator_atom_invariants_(gra))

    # 3. Break ties based on keys.
    idx_dct = relax_class_indices(
        gra, idx_dct, srt_eval_=sort_evaluator_tie_breaking_(gra))

    return idx_dct


def relax_class_indices(gra, idx_dct, srt_eval_):
    """ Relax the class indices for this graph based on some sort value.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param idx_dct: A dictionary mapping atom keys to class indices.
        :type idx_dct: dict
        :param srt_eval_: An evaluator for sort values, based on current class
            indices. Curried such that srt_val_(idx_dct)(key) returns the sort
            value.
    """
    ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    cla_dct = class_dict_from_index_dict(idx_dct)

    # 2. Set up the new_clas list, containing class indices and class keys that
    # are up for re-evaluation.
    new_cla_dct = dict_.filter_by_value(cla_dct, lambda v: len(v) > 1)
    new_clas = sorted(new_cla_dct.items(), reverse=True)

    while new_clas:
        # Pop the next class for re-evaluation
        idx, cla = new_clas.pop(0)

        # Sort and partition the class based on atom invariants and class
        # indices of neighboring atoms.  After the first iteration, the atom
        # invariants have no effect -- only the neighboring class indices cause
        # further subdivision.
        srt_val_ = srt_eval_(idx_dct)
        cla = sorted(cla, key=srt_val_)
        parts = [tuple(ks) for _, ks in itertools.groupby(cla, srt_val_)]

        # Assign new class indices to the partitions and update idx_dct and
        # cla_dct.
        new_idx = idx
        if len(parts) > 1:
            new_idxs = []
            for part in parts:
                cla_dct[new_idx] = part
                idx_dct.update({k: new_idx for k in part})

                # Track the set of nex indices
                new_idxs.append(new_idx)

                # Increment the index for the next class by the number of
                # members in this one, so that that class indices are stable.
                new_idx += len(part)

            # Identify indices of classes with neighboring atoms, as these may
            # be affected by the re-classification.
            ngb_idxs = frozenset()
            for new_idx in new_idxs:
                new_cla = cla_dct[new_idx]

                # Get neighboring keys to the members of this new class.
                ngb_keys = frozenset.union(
                    *map(ngb_keys_dct.__getitem__, new_cla))

                # Get class indices of these neighboring atoms.
                ngb_idxs |= frozenset(map(idx_dct.__getitem__, ngb_keys))

                # Don't revert back to this class, and don't include classes
                # that were already up for re-evaluation.
                ngb_idxs -= {new_idx}
                ngb_idxs -= frozenset(dict(new_clas))

            for ngb_idx in sorted(ngb_idxs):
                ngb_cla = cla_dct[ngb_idx]
                if len(ngb_cla) > 1:
                    new_clas.insert(0, (ngb_idx, cla_dct[ngb_idx]))

    return idx_dct


# Sort evaluators
def sort_evaluator_atom_invariants_(gra):
    """ A sort function based on atom invariants with two levels of currying.

        To get the sort value for a specific key, use
            srt_val = sort_evaluator_atom_invariants_(gra)(idx_dct)(key)

        My reasoning for doing things this way is that `gra` never changes, but
        `idx_dct` does, so we need to be able to update the function with new
        index dictionaries. Ultimately, it is convenient to return a function
        of `key` only because this can be passed to standard python sorting and
        grouping functions such as `sorted()` and `itertools.groupby()`.

        :param gra: molecular graph
        :type gra: automol graph data structure
    """

    def _replace_none(val):
        return numpy.inf if val is None else val

    ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    bnds_dct = atoms_bond_keys(gra)

    anum_dct = atomic_numbers(gra)
    mnum_dct = mass_numbers(gra)
    hnum_dct = atom_implicit_hydrogen_valences(gra)
    apar_dct = dict_.transform_values(atom_stereo_parities(gra), _replace_none)

    bnd_par_dct = bond_stereo_parities(gra)

    def _sortable_bond_stereo_values(bnd_keys):
        bpars = dict_.values_by_key(bnd_par_dct, bnd_keys)
        bpars = tuple(sorted(map(_replace_none, bpars)))
        return bpars

    bpars_dct = dict_.transform_values(bnds_dct, _sortable_bond_stereo_values)

    def _evaluator(idx_dct):
        """ Sort value evaluator based on current class indices.

            :param idx_dct: A dictionary mapping atom keys to class indices.
            :type idx_dct: dict
        """

        def _value(key):
            deg = len(bnds_dct[key])
            anum = anum_dct[key]
            mnum = mnum_dct[key]
            hnum = hnum_dct[key]
            apar = apar_dct[key]
            bpars = bpars_dct[key]
            ngb_idxs = tuple(
                sorted(map(idx_dct.__getitem__, ngb_keys_dct[key])))
            return (deg, anum, mnum, hnum, apar, bpars, ngb_idxs)

        return _value

    return _evaluator


def sort_evaluator_tie_breaking_(gra):
    """ A sort function for tie-breaking with two levels of currying.

        This function is to be called last, after the only remaining classes
        with multiple members are indeed perfectly interchangeable.

        To get the sort value for a specific key, use
            srt_val = sort_evaluator_atom_invariants_(gra)(idx_dct)(key)

        My reasoning for doing things this way is that `gra` never changes, but
        `idx_dct` does, so we need to be able to update the function with new
        index dictionaries. Ultimately, it is convenient to return a function
        of `key` only because this can be passed to standard python sorting and
        grouping functions such as `sorted()` and `itertools.groupby()`.

        :param gra: molecular graph
        :type gra: automol graph data structure
    """

    ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    def _evaluator(idx_dct):
        """ Sort value evaluator based on current class indices.

            :param idx_dct: A dictionary mapping atom keys to class indices.
            :type idx_dct: dict
        """

        def _value(key):
            ngb_idxs = tuple(
                sorted(map(idx_dct.__getitem__, ngb_keys_dct[key])))
            return (key, ngb_idxs)

        return _value

    return _evaluator


# Helpers
def class_dict_from_index_dict(idx_dct):
    """ Obtain a class dictionary from a class index dictionary.

        :param idx_dct: A dictionary mapping atom keys to class indices.
        :type idx_dct: dict
        :returns: A dictionary mapping class indices onto the full set of keys
            for that class, as a sorted tuple.
        :rtype: dict[int: tuple]
    """
    clas = sorted(idx_dct.keys(), key=idx_dct.__getitem__)
    cla_dct = {i: tuple(c)
               for i, c in itertools.groupby(clas, key=idx_dct.__getitem__)}
    return cla_dct


if __name__ == '__main__':
    import automol
#     GEO_STR = """
# C    2.440274   0.058933  -0.692740
# C    1.618916  -0.527269   0.445319
# C    0.152639  -0.210756   0.299288
# C   -0.692099  -1.069400  -0.416694
# C   -2.046782  -0.767410  -0.565130
# C   -2.568232   0.396733  -0.003897
# C   -1.735380   1.260916   0.704387
# C   -0.380457   0.960485   0.853801
# H    2.346592   1.149592  -0.729644
# H    3.498847  -0.186225  -0.560715
# H    2.118103  -0.338675  -1.661135
# H    1.985880  -0.141512   1.404270
# H    1.759752  -1.614397   0.482401
# H   -3.622800   0.631472  -0.120266
# H   -2.140445   2.170539   1.139680
# H   -2.694807  -1.440309  -1.120329
# H    0.258316   1.646462   1.405337
# H   -0.298317  -1.979180  -0.863934
# """
    GEO_STR = """
C     -2.651015   -0.803730    0.673174
C     -1.413459   -0.639646   -0.207699
C     -1.308553    0.802979   -0.732125
C      0.100737    1.208965   -1.205400
C      0.766611    2.166082   -0.212732
C      0.984615   -0.005960   -1.515734
C      1.078078   -1.049159   -0.387249
C      2.374484   -0.923565    0.413812
C     -0.143650   -1.041648    0.544259
H     -3.557123   -0.525575    0.124704
H     -2.589066   -0.176029    1.568637
H     -2.761236   -1.843940    0.997280
H      0.193800    3.097169   -0.139006
H      0.831674    1.741898    0.792568
H      1.778295    2.426720   -0.540834
H      2.438998   -1.714246    1.169043
H      2.444292    0.038985    0.928421
H      3.246275   -1.019060   -0.242139
H     -1.664837    1.512279    0.026061
H     -2.005562    0.903671   -1.574988
H     -0.261382   -2.041315    0.981612
H      0.014703   -0.358487    1.388182
H      0.551474   -0.505544   -2.394009
H      1.983942    0.318361   -1.831221
H     -1.531490   -1.312428   -1.068114
H     -0.014962    1.776522   -2.138918
H      1.114357   -2.033299   -0.875520
"""
    GEO = automol.geom.from_string(GEO_STR)
    print(automol.geom.string(GEO))
    GRA = automol.geom.graph(GEO, stereo=True)
    print(automol.graph.string(GRA))
    GRA = automol.graph.implicit(GRA)
    print(automol.graph.string(GRA))
    # IDX_DCT = canonical_keys(GRA)
    # print("Canonical keys:", IDX_DCT.values())
