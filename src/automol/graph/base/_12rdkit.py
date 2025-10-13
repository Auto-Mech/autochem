"""Interface to RDKit-specific graph functionality."""

import numbers

import rdkit

from ... import util
from ._00core import (
    atom_keys,
    atom_stereo_parities,
    atom_symbols,
    atom_unpaired_electrons,
    bond_keys,
    bond_orders,
    bond_stereo_parities,
    explicit,
    from_data,
    has_stereo,
    set_atom_stereo_parities,
    set_bond_stereo_parities,
    without_bonds_by_orders,
    without_stereo,
)
from ._05stereo import stereocenter_candidates
from ._08canon import (
    bad_stereo_bond_keys_from_kekule,
    from_local_stereo,
    smiles_graph_from_kekule,
    to_local_stereo,
)

ATOM_STEREO_TAG_FROM_BOOL = {
    None: rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
    False: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
    True: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
}
ATOM_STEREO_BOOL_FROM_TAG = dict(map(reversed, ATOM_STEREO_TAG_FROM_BOOL.items()))

BOND_STEREO_TAG_FROM_BOOL = {
    None: rdkit.Chem.rdchem.BondStereo.STEREONONE,
    False: rdkit.Chem.rdchem.BondStereo.STEREOZ,
    True: rdkit.Chem.rdchem.BondStereo.STEREOE,
}
BOND_STEREO_BOOL_FROM_TAG = dict(map(reversed, BOND_STEREO_TAG_FROM_BOOL.items()))
BOND_ORDER_DCT = {
    1: rdkit.Chem.BondType.SINGLE,
    1.5: rdkit.Chem.BondType.ONEANDAHALF,
    2: rdkit.Chem.BondType.DOUBLE,
    3: rdkit.Chem.BondType.TRIPLE,
}
BOND_TYPE_DCT = dict(map(reversed, BOND_ORDER_DCT.items()))


def from_kekule_graph(
    kgr, stereo=False, exp=False, local_stereo=False, label=False, label_dct=None
):
    """Generate an RDKit rdmecule object from a connected molecular graph.

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param exp: Include explicit hydrogens that aren't needed for stereochemistry?
    :type exp: bool
    :param local_stereo: Does the graph have local stereo assignments? defaults to True
    :type local_stereo: bool, optional
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    """
    kgr = without_bonds_by_orders(kgr, ords=[0], skip_dummies=False)
    kgr = smiles_graph_from_kekule(kgr, res_stereo=True, exp=exp, dummy=True)
    rdm, idx_from_key = _from_graph_without_stereo(
        kgr, label=label, label_dct=label_dct
    )

    # If there's not stereo, return early
    if not stereo or not has_stereo(kgr):
        return rdm

    # Otherwise, handle stereo
    kgr = kgr if local_stereo else to_local_stereo(kgr)
    egra = explicit(kgr)
    bad_bkeys = bad_stereo_bond_keys_from_kekule(egra)
    egra = without_stereo(egra, bnd_keys=bad_bkeys)
    erdm, idx_from_key = _from_graph_without_stereo(egra)
    key_from_idx = dict(map(reversed, idx_from_key.items()))
    nkeys_dct = stereocenter_candidates(egra, strict=False)

    # Set atom stereo
    atm_ste_dct = atom_stereo_parities(egra)
    for rda in rdm.GetAtoms():
        idx = rda.GetIdx()
        key = key_from_idx[idx]
        par0 = atm_ste_dct[key]
        if par0 is not None:
            erda = erdm.GetAtoms()[idx]

            nkeys0 = nkeys_dct[key]
            nidxs1 = [b.GetOtherAtomIdx(erda.GetIdx()) for b in erda.GetBonds()]
            nkeys1 = list(map(key_from_idx.__getitem__, nidxs1))

            par1 = util.is_odd_permutation(nkeys0, nkeys1) ^ par0
            rda.SetChiralTag(ATOM_STEREO_TAG_FROM_BOOL[par1])

    # Set bond stereo
    bnd_ste_dct = bond_stereo_parities(kgr)
    for rdb in rdm.GetBonds():
        idxs = (rdb.GetBeginAtomIdx(), rdb.GetEndAtomIdx())
        keys = tuple(map(key_from_idx.get, idxs))
        bkey = frozenset(keys)
        par = bnd_ste_dct[bkey]
        if par is not None:
            bnkeys = nkeys_dct[bkey]
            nkeys1, nkeys2 = [bnkeys[keys.index(k)] for k in sorted(keys)]

            nidx1 = idx_from_key[nkeys1[-1]]
            nidx2 = idx_from_key[nkeys2[-1]]

            rdb.SetStereo(BOND_STEREO_TAG_FROM_BOOL[par])
            rdb.SetStereoAtoms(nidx1, nidx2)

    return rdm


def _from_graph_without_stereo(gra, label=False, label_dct=None):
    """Generate an RDKit rdmecule object from a connected rdmecular graph.

    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    """
    if label_dct is not None:
        label = True

    # Add atoms
    keys = sorted(atom_keys(gra))
    symb_dct = atom_symbols(gra, dummy_symbol="He")
    rad_dct = atom_unpaired_electrons(gra, bond_order=True)

    erdm = rdkit.Chem.EditableMol(rdkit.Chem.Mol())
    for key in keys:
        rda = rdkit.Chem.Atom(symb_dct[key])
        rda.SetNumRadicalElectrons(rad_dct[key])
        erdm.AddAtom(rda)

    # Get key <=> index mappings
    idx_from_key = dict(map(reversed, enumerate(keys)))

    # Add bonds
    bkeys = sorted(bond_keys(gra), key=sorted)
    ord_dct = bond_orders(gra)
    for bkey in bkeys:
        idx1, idx2 = map(idx_from_key.__getitem__, bkey)
        erdm.AddBond(idx1, idx2, BOND_ORDER_DCT[ord_dct[bkey]])

    rdm = erdm.GetMol()

    if label:
        label_dct = {k: k for k in keys} if label_dct is None else label_dct
        # Re-index the label dict to use indices
        label_dct = util.dict_.transform_keys(label_dct, idx_from_key.get)
        for idx, rda in enumerate(rdm.GetAtoms()):
            if idx in label_dct and isinstance(label_dct[idx], numbers.Integral):
                rda.SetProp("molAtomMapNumber", str(label_dct[idx]))

    rdm.UpdatePropertyCache()
    rdkit.Chem.SanitizeMol(rdm)
    return rdm, idx_from_key


def to_graph(exp_rdm, stereo=True, order=False):
    """Generate a connectivity graph from an RDKit molecule object.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :param stereo: Include stereochemistry information?, defaults to True
    :type stereo: bool, optional
    :param order: Include bond orders?, defaults to False
    :type order: bool, optional
    :rtype: automol molecular graph object
    """
    gra = _to_graph_without_stereo(exp_rdm, order=order)

    if not stereo:
        return gra

    exp_rdm = rdkit.Chem.AddHs(exp_rdm)
    exp_gra = _to_graph_without_stereo(exp_rdm)
    nkeys_dct = stereocenter_candidates(exp_gra, strict=False)

    # Assign atom stereo
    for exp_rda in exp_rdm.GetAtoms():
        par0 = ATOM_STEREO_BOOL_FROM_TAG[exp_rda.GetChiralTag()]
        if par0 is not None:
            key = exp_rda.GetIdx()

            nkeys0 = nkeys_dct[key]
            nkeys1 = [b.GetOtherAtomIdx(exp_rda.GetIdx()) for b in exp_rda.GetBonds()]

            par1 = util.is_odd_permutation(nkeys0, nkeys1) ^ par0
            gra = set_atom_stereo_parities(gra, {key: par1})

    # Assign bond stereo
    for exp_rdb in exp_rdm.GetBonds():
        keys = (exp_rdb.GetBeginAtomIdx(), exp_rdb.GetEndAtomIdx())
        bkey = frozenset(keys)

        par0 = BOND_STEREO_BOOL_FROM_TAG[exp_rdb.GetStereo()]
        if par0 is not None:
            bnkeys = nkeys_dct[bkey]
            nkeys1, nkeys2 = [bnkeys[keys.index(k)] for k in sorted(keys)]

            nkey1 = nkeys1[-1]
            nkey2 = nkeys2[-1]
            nkey1_, nkey2_ = exp_rdb.GetStereoAtoms()

            par1 = (nkey1 != nkey1_) ^ (nkey2 != nkey2_) ^ par0
            gra = set_bond_stereo_parities(gra, {bkey: par1})

    gra = from_local_stereo(gra)
    return gra


def _to_graph_without_stereo(rdm, order=False):
    """Generate a connectivity graph from an RDKit molecule object.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :param order: Include bond orders?, defaults to False
    :type order: bool, optional
    :rtype: automol molecular graph object
    """

    def _get_order(rdb):
        return BOND_TYPE_DCT[rdb.GetBondType()] if order else 1

    rdm.UpdatePropertyCache()
    rdkit.Chem.SanitizeMol(rdm)

    rdas = rdm.GetAtoms()
    rdbs = rdm.GetBonds()
    sym_dct = {a.GetIdx(): a.GetSymbol() for a in rdas}
    hyd_dct = {b.GetIdx(): b.GetImplicitValence() for b in rdas}
    ord_dct = {(b.GetBeginAtomIdx(), b.GetEndAtomIdx()): _get_order(b) for b in rdbs}
    gra = from_data(
        atm_symb_dct=sym_dct,
        bnd_keys=ord_dct.keys(),
        atm_imp_hyd_dct=hyd_dct,
        bnd_ord_dct=ord_dct,
    )
    return gra
