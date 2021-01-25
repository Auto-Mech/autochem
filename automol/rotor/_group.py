""" drivers for coordinate scans
"""

import numpy
import automol


# Functions to handle setting up groups and axes used to define torstions
def set_tors_def_info(zma, tors_name, tors_sym, pot,
                      frm_bnd_keys, brk_bnd_keys,
                      rxn_class, saddle=False):
    """ set stuff
    """

    # Set the ts bnd to the break or form keys based on rxn class
    if frm_bnd_keys:
        ts_bnd = frm_bnd_keys
    else:
        ts_bnd = brk_bnd_keys

    group, axis, atm_key = _set_groups_ini(
        zma, tors_name, ts_bnd, saddle)
    if saddle:
        group, axis, pot, chkd_sym_num = _check_saddle_groups(
            zma, rxn_class, group, axis,
            pot, ts_bnd, tors_sym)
    else:
        chkd_sym_num = tors_sym
    group = list(numpy.add(group, 1))
    axis = list(numpy.add(axis, 1))
    if (atm_key+1) != axis[1]:
        axis.reverse()

    return group, axis, pot, chkd_sym_num
