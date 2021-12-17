""" canonicalization functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!

Reference:
Schneider, Sayle, Landrum. J. Chem. Inf. Model. 2015, 55, 10, 2111–2120
"""
import itertools
import functools
import numpy
from phydat import ptab
from automol.util import dict_
import automol.geom.base
from automol.graph.base._core import atom_keys
from automol.graph.base._core import backbone_keys
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import set_atom_stereo_parities
from automol.graph.base._core import set_bond_stereo_parities
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import mass_numbers
from automol.graph.base._core import tetrahedral_atom_keys
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atoms_bond_keys
from automol.graph.base._core import implicit
from automol.graph.base._core import explicit
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import atom_explicit_hydrogen_keys
from automol.graph.base._core import explicit_hydrogen_keys
from automol.graph.base._core import relabel
from automol.graph.base._core import without_bond_orders
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._algo import is_connected
from automol.graph.base._algo import rings_bond_keys
from automol.graph.base._resonance import sp2_bond_keys


# Canonicalization functions
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

        Stereo parities in the graph are assumed to be canonical.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param backbone_only: Consider backbone atoms only?
        :type backbone_only: bool
        :param break_ties: Break ties after keys have been relaxed?
        :type break_ties: bool
        :returns: a dictionary of canonical keys by atom key
        :rtype: dict[int: int]
    """
    can_key_dct, atm_par_dct, bnd_par_dct = canonical_keys_and_stereo_parities(
        gra, backbone_only=backbone_only)

    assert atm_par_dct == atom_stereo_parities(gra), (
        f"Atom stereo parities don't match input. Something is wrong:\n"
        f"input: {atom_stereo_parities(gra)}\n"
        f"return: {atm_par_dct}\n"
    )

    assert bnd_par_dct == bond_stereo_parities(gra), (
        f"Bond stereo parities don't match input. Something is wrong:\n"
        f"input: {bond_stereo_parities(gra)}\n"
        f"return: {bnd_par_dct}\n"
    )

    return can_key_dct


def canonical_keys_and_stereo_parities(gra, backbone_only=True,
                                       atm_par_eval_=None,
                                       bnd_par_eval_=None):
    """ Determine canonical keys and stereo parities for this graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param backbone_only: Consider backbone atoms only?
        :type backbone_only: bool
        :param break_ties: Break ties after keys have been relaxed?
        :type break_ties: bool
        :param atm_par_eval_: An evaluator for atom stereo parities, based on
            the current class indices. Curried such that
            atm_par_val_(idx_dct)(key) returns the sort value.
            If None, the graph is assumed to contain canonical stereo parities,
            and those will be used.
        :param bnd_par_eval_: An evaluator for bond stereo parities, based on
            the current class indices. Curried such that
            bnd_par_val_(idx_dct)(key) returns the sort value.
            If None, the graph is assumed to contain canonical stereo parities,
            and those will be used.
        :returns: a dictionary of canonical keys by atom key
        :rtype: dict[int: int]
    """
    assert is_connected(gra), "Cannot canonicalize disconnected graph."
    assert gra == without_dummy_atoms(gra), (
        "Cannont canonicalize graph with dummy atoms.")

    # By default, assume a graph with canonical stereo.
    atm_par_eval_ = (atom_parity_evaluator_from_canonical_stereo_graph_(gra)
                     if atm_par_eval_ is None else atm_par_eval_)
    bnd_par_eval_ = (bond_parity_evaluator_from_canonical_stereo_graph_(gra)
                     if bnd_par_eval_ is None else bnd_par_eval_)

    # Now, remove stereo parities for the canonicalization
    gra = without_stereo_parities(gra)

    # 1. Calculate initial class indices
    idx_dct = class_indices(gra, break_ties=False)

    # 2. Iterately determine stereo parities.
    old_idx_dct = None
    while idx_dct != old_idx_dct:
        # a. Find stereogenic keys based on the current indices
        atm_keys = list(stereogenic_atom_keys(gra, idx_dct=idx_dct))
        bnd_keys = list(stereogenic_bond_keys(gra, idx_dct=idx_dct))

        # b. Create parity evaluators based on the current indices
        atm_par_ = atm_par_eval_(idx_dct)
        bnd_par_ = bnd_par_eval_(idx_dct)

        # c. Update parities based on the current indices
        atm_par_dct = dict(zip(atm_keys, map(atm_par_, atm_keys)))
        bnd_par_dct = dict(zip(bnd_keys, map(bnd_par_, bnd_keys)))
        gra = set_atom_stereo_parities(gra, atm_par_dct)
        gra = set_bond_stereo_parities(gra, bnd_par_dct)

        # d. Calculate new class indices and return to step a, unless
        #    converged.
        old_idx_dct = idx_dct
        idx_dct = class_indices(gra, idx_dct=idx_dct, break_ties=False)

    # 3. Break ties to determine canonical indices.
    can_key_dct = class_indices(gra, idx_dct=idx_dct,
                                backbone_only=backbone_only, break_ties=True)

    # 4. Read out the final stereo parities.
    atm_par_dct = atom_stereo_parities(gra)
    bnd_par_dct = bond_stereo_parities(gra)

    return can_key_dct, atm_par_dct, bnd_par_dct


def class_indices(gra, backbone_only=True, break_ties=False, idx_dct=None):
    """ Determine symmetry class indices for this graph.

        Unless the `break_ties` flag is turned on, atoms which are
        indistinguishable due to symmetry will end up in the same symmetry
        class.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param backbone_only: Consider backbone atoms only?
        :type backbone_only: bool
        :param break_ties: Break ties after keys have been relaxed?
        :type break_ties: bool
        :param idx_dct: Optionally, pass in initial class indices by key, which
            will be refined.
        :type idx_dct: dict[int: int]
        :returns: a dictionary of class indices by atom key
        :rtype: dict[int: int]
    """
    orig_gra = gra

    # Work with the implicit graph to determine class indices for backbone
    # atoms
    gra = implicit(gra)

    # 1. Initial partition: all atoms in one class, with index 0.
    # idx_dct maps individual atom keys onto their class indices
    idx_dct = {} if idx_dct is None else idx_dct
    # Restrict idx_dct to the backbone keys.
    idx_dct = dict_.by_key(idx_dct, atom_keys(gra), fill_val=0)

    # 2. Relax class indices based on atom invariants
    idx_dct = relax_class_indices(
        gra, idx_dct, srt_eval_=sort_evaluator_atom_invariants_(gra))

    # 3. If requested, break ties based on keys.
    if break_ties:
        idx_dct = relax_class_indices(
            gra, idx_dct, srt_eval_=sort_evaluator_tie_breaking_(gra))

    if not backbone_only:
        idx_dct = add_hydrogen_keys_to_index_dict(orig_gra, idx_dct,
                                                  break_ties=break_ties)

    return idx_dct


# Stereo functions
def stereogenic_atom_keys(gra, idx_dct=None, assigned=False):
    """ Find stereogenic atoms in this graph.

        If the `assigned` flag is set to `False`, only  unassigned stereogenic
        atoms will be detected.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param assigned: Include atoms that already have stereo assignments?
        :param assigned: bool
        :returns: the stereogenic atom keys
        :rtype: frozenset
    """
    # Don't recalculate symmetry classes unless we have to
    idx_dct = class_indices(gra) if idx_dct is None else idx_dct

    gra = without_bond_orders(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in
    idx_dct = add_hydrogen_keys_to_index_dict(gra, idx_dct, break_ties=False)

    atm_keys = tetrahedral_atom_keys(gra)
    if not assigned:
        # Remove assigned stereo keys
        atm_keys -= atom_stereo_keys(gra)

    nkeys_dct = atoms_neighbor_atom_keys(gra)

    def _is_stereogenic(key):
        nkeys = list(nkeys_dct[key])
        idxs = list(map(idx_dct.__getitem__, nkeys))
        return len(set(idxs)) == len(idxs)

    ste_gen_atm_keys = frozenset(filter(_is_stereogenic, atm_keys))
    return ste_gen_atm_keys


def stereogenic_bond_keys(gra, idx_dct=None, assigned=False):
    """ Find stereogenic bonds in this graph.

        If the `assigned` flag is set to `False`, only  unassigned stereogenic
        bonds will be detected.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param assigned: Include bonds that already have stereo assignments?
        :param assigned: bool
        :returns: the stereogenic bond keys
        :rtype: frozenset
    """
    # Don't recalculate symmetry classes unless we have to
    idx_dct = class_indices(gra) if idx_dct is None else idx_dct

    gra = without_bond_orders(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in
    idx_dct = add_hydrogen_keys_to_index_dict(gra, idx_dct, break_ties=False)

    bnd_keys = sp2_bond_keys(gra)
    if not assigned:
        # Remove assigned stereo keys
        bnd_keys -= bond_stereo_keys(gra)

    bnd_keys -= functools.reduce(  # remove double bonds in small rings
        frozenset.union,
        filter(lambda x: len(x) < 8, rings_bond_keys(gra)), frozenset())

    nkeys_dct = atoms_neighbor_atom_keys(gra)

    def _is_stereogenic(key):

        def _is_asymmetric_on_bond(atm1_key, atm2_key):
            nkeys = list(nkeys_dct[atm1_key] - {atm2_key})

            if not nkeys:                # C=:O:
                # Atoms without neighbors are automatically symmetric
                ret = False
            elif len(nkeys) == 1:        # C=N:-X
                # Atoms without 1 neighbor are automatically asymmetric
                ret = True
            else:
                # For atoms with 2 neighbors, we need to determine whether or
                # not they are symmetric from the class indices.
                assert len(nkeys) == 2   # C=C(-X)-Y
                ret = idx_dct[nkeys[0]] != idx_dct[nkeys[1]]

            return ret

        atm1_key, atm2_key = key
        return (_is_asymmetric_on_bond(atm1_key, atm2_key) and
                _is_asymmetric_on_bond(atm2_key, atm1_key))

    ste_gen_bnd_keys = frozenset(filter(_is_stereogenic, bnd_keys))
    return ste_gen_bnd_keys


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


# Parity evaluators
def atom_parity_evaluator_from_geometry_(gra, geo, geo_idx_dct=None):
    r""" A stereo parity evaluator for atoms in a geometry

        Parity is defined as follows:

        The four keys passed in are apices of a tetrahedron. Looking at 2, 3,
        and 4 from 1, they will either ascend in clockwise or counterclockwise
        order.

        If ascending in clockwise order, the parity is True ('+').
        If ascending in counterclockwise order, the parity is False ('-').

             2              2
            /1\            /1\
           4---3          3---4
           clockwise      counterclockwise
           True           False
           '+'            '-'

        (Viewed looking down from 1)

        If only three keys are passed in, they will be treated as keys 2, 3,
        and 4 above and it will be assumed that there is a lone pair at 1.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
        :returns: A parity evaluator, curried such that par_eval_(idx_dct)(key)
            returns the parity for a given atom, given a set of class indices.
    """
    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(atm_keys)})

    xyzs = automol.geom.base.coordinates(geo)
    xyz_dct = {k: xyzs[geo_idx_dct[k]] for k in atm_keys}

    def _evaluator(idx_dct):
        """ Parity evaluator based on current class indices.

            :param idx_dct: A dictionary mapping atom keys to class indices.
            :type idx_dct: dict
        """
        idx_dct = add_hydrogen_keys_to_index_dict(gra, idx_dct,
                                                  break_ties=False)

        def _parity(key, keys):
            # Sort the keys by class index
            keys = sorted(keys, key=idx_dct.__getitem__)

            # If there are only three groups, use the stereo atom itself as the
            # top apex of the tetrahedron.
            if len(keys) != 4:
                assert len(keys) == 3
                keys = [key] + list(keys)

            xyzs = list(map(list, map(xyz_dct.__getitem__, keys)))
            det_mat = numpy.ones((4, 4))
            det_mat[:, 1:] = xyzs
            det_val = numpy.linalg.det(det_mat)
            assert det_val != 0.  # for now, assume no four-atom planes
            par = det_val > 0.
            return par

        return _parity

    return _evaluator


def atom_parity_evaluator_from_canonical_stereo_graph_(gra):
    """ A stereo parity evaluator atoms in a graph with canonical stereo

        Stereo is already assumed to be canonical, so the class indices passed
        in are ignored.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: A parity evaluator, curried such that par_eval_(idx_dct)(key)
            returns the parity for a given atom, given a set of class indices.
    """
    atm_par_dct = atom_stereo_parities(gra)

    def _evaluator(idx_dct):
        """ Parity evaluator based on current class indices

            Class indices are ignored, since the parities are assumed to be
            canonical.

            :param idx_dct: A dictionary mapping atom keys to class indices.
            :type idx_dct: dict
        """
        # Do-nothing line to prevent linting complaint
        assert idx_dct or not idx_dct

        def _parity(key):
            return atm_par_dct[key]

        return _parity

    return _evaluator


def bond_parity_evaluator_from_canonical_stereo_graph_(gra):
    """ A stereo parity evaluator bonds in a graph with canonical stereo

        Stereo is already assumed to be canonical, so the class indices passed
        in are ignored.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: A parity evaluator, curried such that par_eval_(idx_dct)(key)
            returns the parity for a given bond, given a set of class indices.
    """
    bnd_par_dct = bond_stereo_parities(gra)

    def _evaluator(idx_dct):
        """ Parity evaluator based on current class indices

            Class indices are ignored, since the parities are assumed to be
            canonical.

            :param idx_dct: A dictionary mapping bond keys to class indices.
            :type idx_dct: dict
        """
        # Do-nothing line to prevent linting complaint
        assert idx_dct or not idx_dct

        def _parity(key):
            return bnd_par_dct[key]

        return _parity

    return _evaluator


# Canonicalization helpers
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


def add_hydrogen_keys_to_index_dict(gra, idx_dct, break_ties=False):
    """ Add explicit hydrogen keys to the class index dictionary
    """
    idx_dct = idx_dct.copy()

    hyd_keys_pool = explicit_hydrogen_keys(gra)

    if hyd_keys_pool:
        bbn_keys = sorted(backbone_keys(gra), key=idx_dct.__getitem__)
        gra = explicit(gra)

        hyd_keys_dct = atom_explicit_hydrogen_keys(gra)

        next_idx = len(bbn_keys)

        # If not breaking ties, partition hydrogens bonded to atoms from
        # the same class into classes and use corresponding indices.
        hyd_idx_dct = {}
        if not break_ties:
            # Partition hydrogens into classes based on their parent atoms.
            bbn_parts = sorted_classes_from_index_dict(idx_dct)
            get_hyd_key_ = hyd_keys_dct.__getitem__
            hyd_parts = [
                list(itertools.chain(*map(get_hyd_key_, bbn_part)))
                for bbn_part in bbn_parts]
            for hyd_part in hyd_parts:
                hyd_idx_dct.update({k: next_idx for k in hyd_part})
                next_idx += len(hyd_part)
        # Otherwise, give each hydrogen a unique label.
        else:
            srt_hyd_keys = itertools.chain(
                *map(hyd_keys_dct.__getitem__, bbn_keys))
            hyd_idx_dct.update({k: i+next_idx for i, k in
                                enumerate(srt_hyd_keys)})

        hyd_idx_dct = dict_.by_key(hyd_idx_dct, hyd_keys_pool)
        idx_dct.update(hyd_idx_dct)

    return idx_dct


def class_dict_from_index_dict(idx_dct):
    """ Obtain a class dictionary from a class index dictionary.

        :param idx_dct: A dictionary mapping atom keys to class indices.
        :type idx_dct: dict
        :returns: A dictionary mapping class indices onto the full set of keys
            for that class, as a sorted tuple.
        :rtype: dict[int: tuple]
    """
    keys = sorted(idx_dct.keys())
    clas = sorted(keys, key=idx_dct.__getitem__)
    cla_dct = {i: tuple(c)
               for i, c in itertools.groupby(clas, key=idx_dct.__getitem__)}
    return cla_dct


def sorted_classes_from_index_dict(idx_dct):
    """ Obtain classes from index dict, sorted by class index.

        :param idx_dct: A dictionary mapping atom keys to class indices.
        :type idx_dct: dict
        :returns: A tuple of tuples of keys for each class, sorted by class
            index.
        :rtype: tuple[tuple[int]]
    """
    keys = sorted(idx_dct.keys())
    clas = sorted(keys, key=idx_dct.__getitem__)
    cla_dct = tuple(
        tuple(c) for _, c in itertools.groupby(clas, key=idx_dct.__getitem__))
    return cla_dct


# if __name__ == '__main__':
#     import automol
#     # ICH = automol.smiles.inchi('C[C@](F)(N)O')
#     # ICH = automol.smiles.inchi('N(Br)(Cl)F')
#     # ICH = automol.smiles.inchi(r'[H]/N=C\C(\C=N\[H])=N\[H]')
#     ICH = automol.smiles.inchi('F[C@@H]([C@@H](F)Cl)[C@H](F)Cl')
#     print(ICH)
#     GEO = automol.inchi.geometry(ICH)
#     print(automol.geom.string(GEO))
#     GRA = automol.geom.graph(GEO, stereo=True)
#     canonical_keys_and_stereo_parities(
#         GRA,
#         atm_par_eval_=atom_parity_evaluator_from_canonical_stereo_graph_(GRA),
#         bnd_par_eval_=bond_parity_evaluator_from_canonical_stereo_graph_(GRA),
#     )
