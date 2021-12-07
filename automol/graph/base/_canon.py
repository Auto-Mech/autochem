""" canonicalization functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
import numpy
from automol.util import dict_
from automol.graph.base._core import atom_keys
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import implicit
from automol.graph.base._core import atomic_numbers
from automol.graph.base._core import mass_numbers
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atoms_bond_keys


def canonical_keys(gra):
    """ Determine canonical keys for this graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :rtype: dict[int: int]
    """
    gra = implicit(gra)
    ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    # 1. Initial partition: all atoms in one class, with index 0
    # idx_dct maps individual atom keys onto their class indices
    # cla_dct maps class indices onto the full set of keys for that class,
    # stored as a tuple.
    idx_dct = dict_.by_key({}, atom_keys(gra), fill_val=0)
    cla_dct = {0: tuple(atom_keys(gra))}

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
        srt_val_ = sort_value_atom_invariants_(gra, idx_dct)
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

            print(sorted(ngb_idxs))
            for ngb_idx in sorted(ngb_idxs):
                ngb_cla = cla_dct[ngb_idx]
                if len(ngb_cla) > 1:
                    new_clas.insert(0, (ngb_idx, cla_dct[ngb_idx]))

        print()
        print("IDX_DCT", idx_dct.values())
        print("CLA_DCT", cla_dct)
        print("NEW_CLAS", new_clas)


def canonical_keys_old(gra):
    """ Determine canonical keys for this graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :rtype: dict[int: int]
    """
    # orig_gra = gra

    gra = implicit(gra)
    ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    # 1. INITIAL PARTITION: ALL ATOMS
    idx_dct = dict_.by_key({}, atom_keys(gra), fill_val=0)
    cla_dct = {0: atom_keys(gra)}

    # 1B. SET UP NEW_CLAS (CLASSES UP FOR RE-EVALUATION)
    new_clas = dict_.filter_by_value(cla_dct, lambda v: len(v) > 1).items()
    new_clas = sorted(new_clas, reverse=True)

    # 1C. POP THE FIRST (ONLY) CLASS IN NEW_CLAS
    idx, cla = new_clas.pop(0)

    # 2. GROUP/SORT BY INITIAL ATOM INVARIANT
    srt_val_ = sort_value_atom_invariants_(gra, idx_dct)
    res = sorted(
        (v, list(ks)) for v, ks in itertools.groupby(cla, srt_val_))

    # 2B. UPDATE CLASS INDICES AND RECORD NEIGHBORS
    new_idx = idx
    ngb_keys = frozenset()
    for _, new_cla in res:
        cla_dct[new_idx] = new_cla
        idx_dct.update({k: new_idx for k in new_cla})

        # Set the next index
        new_idx += len(new_cla)

        # Update the list of neighbors potentially affected by this change
        ngb_keys |= frozenset.union(*map(ngb_keys_dct.__getitem__, new_cla))

    print(idx_dct.values())
    print(cla_dct)
    print(ngb_keys)

    # 2C. IDENTIFY CLASSES UP FOR RE-EVALUATION
    ngb_idxs = frozenset(map(idx_dct.__getitem__, ngb_keys))

    # 2D. SET UP NEW_CLAS (CLASSES UP FOR RE-EVALUATION)
    new_cla_dct = dict_.by_key(cla_dct, ngb_idxs)
    new_clas = dict_.filter_by_value(new_cla_dct, lambda v: len(v) > 1).items()
    print(new_clas)
    new_clas = sorted(new_clas, reverse=True)
    print(new_clas)

    old_clas = None

    # 3. RELAXATION: CONTINUE TO UPDATE BASED ON NEIGHBORS
    while new_clas != old_clas:
        print('NEW', new_clas)
        print('OLD', new_clas)
        idx, cla = new_clas.pop(0)
        print(idx, cla)

        # 3A. GROUP/SORT BY NEIGHBOR CLASS INDICES
        srt_val_ = sort_value_atom_invariants_(gra, idx_dct)
        cla = sorted(cla, key=srt_val_)
        grps = [(v, tuple(ks)) for v, ks in itertools.groupby(cla, srt_val_)]

        # 3B. UPDATE CLASS INDICES AND RECORD NEIGHBORS
        new_idx = idx
        ngb_keys = frozenset()
        for _, new_cla in grps:
            cla_dct[new_idx] = new_cla
            idx_dct.update({k: new_idx for k in new_cla})

            # Set the next index
            new_idx += len(new_cla)

            # Update the list of neighbors potentially affected by this change
            ngb_keys |= frozenset.union(
                *map(ngb_keys_dct.__getitem__, new_cla))

        print(idx_dct.values())
        print(cla_dct)
        print(ngb_keys)
        # IDENTIFY CLASSES UP FOR RE-EVALUATION AND SET UP NEW_CLAS
        ngb_idxs = frozenset(map(idx_dct.__getitem__, ngb_keys))
        new_cla_dct = dict_.by_key(cla_dct, ngb_idxs)
        new_cla_dct = dict_.filter_by_value(new_cla_dct, lambda v: len(v) > 1)

        old_clas = new_clas
        new_clas = sorted(new_cla_dct.items(), reverse=True) + new_clas


def sort_value_atom_invariants_(gra, idx_dct):
    """ Returns a function returning sort values based on atom invariants
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

    def _value(key):
        deg = len(bnds_dct[key])
        anum = anum_dct[key]
        mnum = mnum_dct[key]
        hnum = hnum_dct[key]
        apar = apar_dct[key]
        bpars = bpars_dct[key]
        ngb_idxs = tuple(sorted(map(idx_dct.__getitem__, ngb_keys_dct[key])))
        return (deg, anum, mnum, hnum, apar, bpars, ngb_idxs)

    return _value


if __name__ == '__main__':
    import automol
    GEO_STR = """
C    2.440274   0.058933  -0.692740
C    1.618916  -0.527269   0.445319
C    0.152639  -0.210756   0.299288
C   -0.692099  -1.069400  -0.416694
C   -2.046782  -0.767410  -0.565130
C   -2.568232   0.396733  -0.003897
C   -1.735380   1.260916   0.704387
C   -0.380457   0.960485   0.853801
H    2.346592   1.149592  -0.729644
H    3.498847  -0.186225  -0.560715
H    2.118103  -0.338675  -1.661135
H    1.985880  -0.141512   1.404270
H    1.759752  -1.614397   0.482401
H   -3.622800   0.631472  -0.120266
H   -2.140445   2.170539   1.139680
H   -2.694807  -1.440309  -1.120329
H    0.258316   1.646462   1.405337
H   -0.298317  -1.979180  -0.863934
"""
    GEO = automol.geom.from_string(GEO_STR)
    GRA = automol.geom.graph(GEO)
    GRA = automol.graph.implicit(GRA)
    canonical_keys(GRA)
