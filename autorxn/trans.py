""" graph transformations
"""
import itertools
import numbers
import automol


def bond_formation_keys(tra):
    """ bond formation keys
    """
    form_keys, _, _ = tra
    return form_keys


def bond_break_keys(tra):
    """ bond break keys
    """
    _, break_keys, _ = tra
    return break_keys


def product_atom_keys(tra):
    """ atom keys in the product, by atom key
    """
    _, _, atm_key_dct = tra
    return atm_key_dct


def apply(tra, cgr):
    """ apply a transformation to a molecular graph
    """
    assert cgr == automol.graph.without_bond_orders(cgr)
    assert cgr == automol.graph.without_stereo_parities(cgr)

    form_keys, break_keys, iso_dct = tra
    for form_key in form_keys:
        if _is_atom_key(form_key):
            cgr = automol.graph.add_atom_implicit_hydrogen_valences(
                cgr, {form_key: +1})
        else:
            assert _is_bond_key(form_key)
            assert form_key not in automol.graph.bond_keys(cgr)
            cgr = automol.graph.add_bonds(cgr, [form_key])
    for break_key in break_keys:
        if _is_atom_key(break_key):
            cgr = automol.graph.add_atom_implicit_hydrogen_valences(
                cgr, {break_key: -1})
        else:
            assert _is_bond_key(break_key)
            assert break_key in automol.graph.bond_keys(cgr)
            cgr = automol.graph.remove_bonds(cgr, [break_key])
    cgr = automol.graph.relabel(cgr, iso_dct)
    return cgr


def invert(tra):
    """ invert a transformation
    """
    orig_form_keys, orig_break_keys, iso_dct = tra
    form_keys = frozenset(
        frozenset(map(iso_dct.__getitem__, key)) if not _is_atom_key(key)
        else iso_dct[key] for key in orig_break_keys)
    break_keys = frozenset(
        frozenset(map(iso_dct.__getitem__, key)) if not _is_atom_key(key)
        else iso_dct[key] for key in orig_form_keys)
    iso_dct = dict(map(reversed, iso_dct.items()))
    tra = (form_keys, break_keys, iso_dct)
    return tra


def direct_sum(tra1, tra2):
    """ direct sum of two orthogonal transformations
    """
    form_keys1, break_keys1, iso_dct1 = tra1
    form_keys2, break_keys2, iso_dct2 = tra2
    assert not set(form_keys1) & set(form_keys2)
    assert not set(break_keys1) & set(break_keys2)
    assert not set(iso_dct1.keys()) & set(iso_dct2.keys())
    assert not set(iso_dct1.values()) & set(iso_dct2.values())
    form_keys = set(form_keys1) | set(form_keys2)
    break_keys = set(break_keys1) | set(break_keys2)
    iso_dct = {}
    iso_dct.update(iso_dct1)
    iso_dct.update(iso_dct2)
    tra = (form_keys, break_keys, iso_dct)
    return tra


def stereo_compatible_products(sgr, tra):
    """ expand stereo for the products graph of a transformation
    """
    sgr1 = sgr
    assert not automol.graph.stereogenic_atom_keys(sgr1)
    assert not automol.graph.stereogenic_bond_keys(sgr1)

    cgr1 = automol.graph.without_stereo_parities(sgr1)
    cgr2 = apply(tra, cgr1)

    sgrs2 = tuple(filter(lambda x: is_stereo_compatible(tra, sgr1, x),
                         automol.graph.stereomers(cgr2)))
    return sgrs2


def is_stereo_compatible(tra, sgr1, sgr2):
    """ is this transformation compatible with the stereo assignments for these
    molecular graph
    """
    cgr1 = automol.graph.without_stereo_parities(sgr1)
    cgr2 = automol.graph.without_stereo_parities(sgr2)
    assert apply(tra, cgr1) == cgr2

    # use a common set of labels for testing compatibility
    sgr1 = automol.graph.relabel(sgr1, product_atom_keys(tra))
    atm_ste_keys1 = automol.graph.atom_stereo_keys(sgr1)
    bnd_ste_keys1 = automol.graph.bond_stereo_keys(sgr1)
    atm_ste_keys2 = automol.graph.atom_stereo_keys(sgr2)
    bnd_ste_keys2 = automol.graph.bond_stereo_keys(sgr2)

    atm_ste_keys = atm_ste_keys1 & atm_ste_keys2
    bnd_ste_keys = bnd_ste_keys1 & bnd_ste_keys2

    atm_ste_par_dct1 = automol.graph.atom_stereo_parities(sgr1)
    bnd_ste_par_dct1 = automol.graph.bond_stereo_parities(sgr1)
    atm_ste_par_dct2 = automol.graph.atom_stereo_parities(sgr2)
    bnd_ste_par_dct2 = automol.graph.bond_stereo_parities(sgr2)

    atm_ngb_keys_dct1 = automol.graph.atom_neighbor_keys(sgr1)
    atm_ngb_keys_dct2 = automol.graph.atom_neighbor_keys(sgr2)

    ret = True

    for atm_key in atm_ste_keys:
        par1 = atm_ste_par_dct1[atm_key]
        par2 = atm_ste_par_dct2[atm_key]

        atm_ngb_keys1 = automol.graph.stereo_sorted_atom_neighbor_keys(
            sgr1, atm_key, atm_ngb_keys_dct1[atm_key])
        atm_ngb_keys2 = automol.graph.stereo_sorted_atom_neighbor_keys(
            sgr2, atm_key, atm_ngb_keys_dct2[atm_key])

        if _permutation_parity(atm_ngb_keys1, atm_ngb_keys2):
            ret &= (par1 == par2)
        else:
            ret &= (par1 != par2)

    for bnd_key in bnd_ste_keys:
        par1 = bnd_ste_par_dct1[bnd_key]
        par2 = bnd_ste_par_dct2[bnd_key]

        atm1_key, atm2_key = bnd_key

        atm1_ngb_key1 = automol.graph.stereo_sorted_atom_neighbor_keys(
            sgr1, atm1_key, atm_ngb_keys_dct1[atm1_key] - {atm2_key})[0]
        atm2_ngb_key1 = automol.graph.stereo_sorted_atom_neighbor_keys(
            sgr1, atm2_key, atm_ngb_keys_dct1[atm2_key] - {atm1_key})[0]
        atm1_ngb_key2 = automol.graph.stereo_sorted_atom_neighbor_keys(
            sgr2, atm1_key, atm_ngb_keys_dct2[atm1_key] - {atm2_key})[0]
        atm2_ngb_key2 = automol.graph.stereo_sorted_atom_neighbor_keys(
            sgr2, atm2_key, atm_ngb_keys_dct2[atm2_key] - {atm1_key})[0]

        if not ((atm1_ngb_key1 != atm1_ngb_key2) ^
                (atm2_ngb_key1 != atm2_ngb_key2)):
            ret &= (par1 == par2)
        else:
            ret &= (par1 != par2)

    return ret


def migration(gra1, gra2):
    """ migration transformation
    """
    tra = None
    rct_gras = automol.graph.connected_components(gra1)
    prd_gras = automol.graph.connected_components(gra2)

    if len(rct_gras) == 1 and len(prd_gras) == 1:
        r_gra, = rct_gras
        p_gra, = prd_gras

        r_atm_keys = automol.graph.resonance_dominant_radical_atom_keys(r_gra)
        p_atm_keys = automol.graph.resonance_dominant_radical_atom_keys(p_gra)
        for r_atm_key, p_atm_key in itertools.product(r_atm_keys, p_atm_keys):
            rh_gra = automol.graph.add_atom_implicit_hydrogen_valences(
                r_gra, {r_atm_key: +1})
            ph_gra = automol.graph.add_atom_implicit_hydrogen_valences(
                p_gra, {p_atm_key: +1})
            iso_dct = automol.graph.backbone_isomorphism(rh_gra, ph_gra)
            if iso_dct:
                inv_iso_dct = dict(map(reversed, iso_dct.items()))
                form_keys = (r_atm_key,)
                break_keys = (inv_iso_dct[p_atm_key],)
                tra = (form_keys, break_keys, iso_dct)
                break
    return tra


def addition(gra1, gra2):
    """ addition transformation
    """
    tra = None
    rct_gras = automol.graph.connected_components(gra1)
    prd_gras = automol.graph.connected_components(gra2)

    if len(rct_gras) == 2 and len(prd_gras) == 1:
        x_gra, y_gra = rct_gras
        xy_gra, = prd_gras
        assert not (automol.graph.atom_keys(x_gra) &
                    automol.graph.atom_keys(y_gra))
        u_gra = automol.graph.union(x_gra, y_gra)
        x_atm_keys = automol.graph.unsaturated_atom_keys(x_gra)
        y_atm_keys = automol.graph.unsaturated_atom_keys(y_gra)
        for x_atm_key, y_atm_key in itertools.product(x_atm_keys, y_atm_keys):
            xy_gra_ = automol.graph.add_bonds(u_gra, [{x_atm_key, y_atm_key}])
            iso_dct = automol.graph.backbone_isomorphism(xy_gra_, xy_gra)
            if iso_dct:
                form_keys = (frozenset({x_atm_key, y_atm_key}),)
                break_keys = ()
                tra = (form_keys, break_keys, iso_dct)
                break

    return tra


def abstraction(gra1, gra2):
    """ abstraction transformation
    """
    tra = None
    rct_gras = automol.graph.connected_components(gra1)
    prd_gras = automol.graph.connected_components(gra2)

    if len(rct_gras) == len(prd_gras) == 2:
        for r_gras in itertools.permutations(rct_gras):
            for p_gras in itertools.permutations(prd_gras):
                r_fmls = list(map(automol.graph.formula, r_gras))
                p_fmls = list(map(automol.graph.formula, p_gras))
                if (automol.formula.add_hydrogen(p_fmls[0]) == r_fmls[0] and
                        automol.formula.add_hydrogen(r_fmls[1]) == p_fmls[1]):
                    q1h_gra, q2_gra = r_gras
                    q1_gra, q2h_gra = p_gras
                    tra1 = _partial_abstraction(q1_gra, q1h_gra)
                    tra2 = _partial_abstraction(q2_gra, q2h_gra)
                    if tra1 is not None and tra2 is not None:
                        tra = direct_sum(invert(tra1), tra2)

    return tra


def _partial_abstraction(q_gra, qh_gra):
    """ backward abstraction transformation
    """
    tra = None
    atm_keys = automol.graph.resonance_dominant_radical_atom_keys(q_gra)
    for atm_key in atm_keys:
        qh_gra_ = automol.graph.add_atom_implicit_hydrogen_valences(
            q_gra, {atm_key: +1})
        iso_dct = automol.graph.backbone_isomorphism(qh_gra_, qh_gra)
        if iso_dct:
            form_keys = (atm_key,)
            break_keys = ()
            tra = (form_keys, break_keys, iso_dct)
            break

    return tra


def _is_atom_key(key):
    return isinstance(key, numbers.Integral)


def _is_bond_key(key):
    return isinstance(key, frozenset) and len(key) == 2 and all(
        isinstance(item, numbers.Integral) for item in key)


def _permutation_parity(seq, ref_seq):
    size = len(seq)
    assert sorted(seq) == sorted(ref_seq) and len(set(seq)) == size
    perm = [ref_seq.index(val) for val in seq]

    sgn = 1
    for idx in range(size):
        if perm[idx] != idx:
            sgn *= -1
            swap_idx = perm.index(idx)
            perm[idx], perm[swap_idx] = perm[swap_idx], perm[idx]

    par = (sgn == 1)

    return par
