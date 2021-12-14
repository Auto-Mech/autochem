""" canonicalization functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!

Reference:
Schneider, Sayle, Landrum. J. Chem. Inf. Model. 2015, 55, 10, 2111–2120
"""
import itertools
import numpy
from phydat import ptab
from automol.util import dict_
from automol.graph.base._core import atom_keys
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import mass_numbers
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import atom_explicit_hydrogen_keys
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atoms_bond_keys
from automol.graph.base._core import implicit
from automol.graph.base._core import relabel
from automol.graph.base._algo import is_connected


def canonical(gra):
    """ A graph relabeled with canonical keys

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: a new molecular graph with canonical keys; if explicit
            hydrogens are included, they will be relabeled as well
    """
    can_key_dct = canonical_keys(gra, backbone_only=False)
    return relabel(gra, can_key_dct)


def canonical_keys(gra, backbone_only=True):
    """ Determine canonical keys for this graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param backbone_only: Consider backbone atoms only?
        :type backbone_only: bool
        :returns: a dictionary of canonical keys by atom key
        :rtype: dict[int: int]
    """
    assert is_connected(gra), "Cannot canonicalize disconnected graph."

    gra = without_dummy_atoms(gra)
    orig_gra = gra

    # Work with the implicit graph to determine canonical keys for backbone
    # atoms
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

    if not backbone_only:
        bbn_keys = atom_keys(gra)
        offset = len(bbn_keys)

        exp_hyd_keys_dct = atom_explicit_hydrogen_keys(orig_gra)
        srt_exp_hyd_keys = [k
                            for b in sorted(bbn_keys, key=idx_dct.__getitem__)
                            for k in sorted(exp_hyd_keys_dct[b])]
        idx_dct.update({k: i+offset for i, k in enumerate(srt_exp_hyd_keys)})

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

        The canonical order has been modified relative to the one used by
        Scheider (2015) to be more InChI-like (following Hill ordering).

        :param gra: molecular graph
        :type gra: automol graph data structure
    """

    def _replace_none(val):
        return numpy.inf if val is None else val

    def _hill_normalize_symbol(symb):
        """ normalize atomic symbols to make them sort according to the hill
            system
            (C first, H second, others in alphabetical order)
        """
        symb = ptab.to_symbol(symb)
        if symb == 'C':
            symb = ''
        if symb == 'H':
            symb = '1'
        return symb

    ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    bnds_dct = atoms_bond_keys(gra)

    symb_dct = dict_.transform_values(
        atom_symbols(gra), _hill_normalize_symbol)
    hnum_dct = atom_implicit_hydrogen_valences(gra)
    mnum_dct = mass_numbers(gra)
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
            symb = symb_dct[key]        # symbol
            deg = len(bnds_dct[key])    # number of bonds
            hnum = hnum_dct[key]        # number of hydrogens
            mnum = mnum_dct[key]
            apar = apar_dct[key]
            bpars = bpars_dct[key]
            ngb_idxs = tuple(
                sorted(map(idx_dct.__getitem__, ngb_keys_dct[key])))
            return (symb, deg, hnum, mnum, apar, bpars, ngb_idxs)

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
            return (ngb_idxs, key)

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
