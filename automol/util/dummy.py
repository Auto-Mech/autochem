"""
  Useful functions for dealing with dummy atom shifts
"""

import automol.zmat


# Add the downshift function from rotor
def shift_down(zma, vals):
    """
    Build a remdummy list that tells how to shift the groups
    """

    dummy_idxs = sorted(automol.zmat.atom_indices(zma, 'X', match=True))

    if dummy_idxs:
        remdummy = [0 for _ in range(automol.zmat.count(zma))]
        for dummy in dummy_idxs:
            for idx, _ in enumerate(remdummy):
                if dummy < idx:
                    remdummy[idx] += 1

        vals1 = tuple(val+1 for val in vals)
        vals2 = tuple(val-remdummy[val-1] for val in vals1)
        shift_vals = tuple(val-1 for val in vals2)

    else:
        shift_vals = vals

    return shift_vals


def shift_up(zma, idxs):
    """ shift up from the dummy idxs
    """

    dummy_idxs = sorted(automol.zmat.atom_indices(zma, 'X', match=True))

    shift_idxs = []
    for idx in idxs:
        new_idx = idx
        for dummy_idx in dummy_idxs:
            if idx >= dummy_idx:
                new_idx += 1
        shift_idxs.append(new_idx)

    return tuple(shift_idxs)
