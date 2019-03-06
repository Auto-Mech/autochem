""" miscellaneous functions
"""
from ._dict import by_key as _by_key
from ._expl import implicit as _implicit
from ._core import bond_keys as _bond_keys
from ._core import (atom_implicit_hydrogen_valences as
                    _atom_implicit_hydrogen_valences)


def bond_symmetry_numbers(xgr):
    """ symmetry numbers, by bond

    the (approximate) symmetry number of the torsional potential for this bond,
    based on the hydrogen counts for each atom
    """
    imp_xgr = _implicit(xgr)
    atm_imp_hyd_vlc_dct = _atom_implicit_hydrogen_valences(imp_xgr)

    bnd_keys = _bond_keys(imp_xgr)
    bnd_max_hyd_vlcs = [max(map(atm_imp_hyd_vlc_dct.__getitem__, bnd_key))
                        for bnd_key in bnd_keys]
    bnd_sym_nums = [3 if vlc == 3 else 1 for vlc in bnd_max_hyd_vlcs]
    bnd_sym_num_dct = dict(zip(bnd_keys, bnd_sym_nums))

    # fill in the rest of the bonds for completeness
    bnd_sym_num_dct = _by_key(bnd_sym_num_dct, _bond_keys(xgr), fill_val=1)
    return bnd_sym_num_dct
