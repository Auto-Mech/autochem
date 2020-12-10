"""
  Parameters
"""

import inspect
import itertools


class REACTION_CLASS():
    """ Names of supported reaction classes
    """
    TRIVIAL = 'trivial'
    # Unimolecular reactions
    BETA_SCISSION = 'beta scission'
    ELIMINATION = 'elimination'
    RING_FORM_SCISSION = 'ring forming scission'
    # Bimolecular reactions
    HYDROGEN_MIGRATION = 'hydrogen migration'
    HYDROGEN_ABSTRACTION = 'hydrogen abstraction'
    ADDITION = 'addition'
    INSERTION = 'insertion'
    SUBSTITUTION = 'substitution'


def is_reation_class(rclass):
    """ Check if class in list of REACTION CLASS
    """
    return rclass in _values(REACTION_CLASS)


def _values(cls):
    """ list the values of a parameter class
    """
    assert inspect.isclass(cls)
    vals = tuple(val for val in _public_attributes(cls)
                 if not inspect.isclass(val))
    return vals


def all_values(cls):
    """ recursively list the values of a parameter class tree
    """
    assert inspect.isclass(cls)
    vals = tuple(itertools.chain(*(
        [val] if not inspect.isclass(val) else all_values(val)
        for val in _public_attributes(cls))))
    return vals


def _public_attributes(cls):
    return tuple(val for name, val in
                 inspect.getmembers(cls, lambda x: not inspect.isroutine(x))
                 if not name.startswith('_') and not inspect.isfunction(val))
