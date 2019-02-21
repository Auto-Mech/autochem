""" not sure what to call this yet
"""
from functools import partial as _partial


def unique(itms, equiv):
    """ unique items from a list, according to binary comparison `equiv`
    """
    uniq_itms = []
    for itm in itms:
        if not any(map(_partial(equiv, itm), uniq_itms)):
            uniq_itms.append(itm)

    return tuple(uniq_itms)
