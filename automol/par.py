""" common automol parameters
"""
import inspect
import itertools


class ReactionClass:
    """ Names of supported reaction classes
    """
    TRIVIAL = 'trivial'
    # Unimolecular reactions
    HYDROGEN_MIGRATION = 'hydrogen migration'
    BETA_SCISSION = 'beta scission'
    ELIMINATION = 'elimination'
    RING_FORM_SCISSION = 'ring forming scission'
    # Bimolecular reactions
    HYDROGEN_ABSTRACTION = 'hydrogen abstraction'
    ADDITION = 'addition'
    INSERTION = 'insertion'
    SUBSTITUTION = 'substitution'


REVERSE_REACTION_DCT = {
    ReactionClass.HYDROGEN_MIGRATION: ReactionClass.HYDROGEN_MIGRATION,
    ReactionClass.HYDROGEN_ABSTRACTION: ReactionClass.HYDROGEN_ABSTRACTION,
    ReactionClass.ADDITION: ReactionClass.BETA_SCISSION,
    ReactionClass.BETA_SCISSION: ReactionClass.ADDITION,
    # ReactionClass.ELIMINATION: ?
    # ReactionClass.INSERTION: ?
    # ReactionClass.SUBSTITUTION: ?
}


def is_reaction_class(rxn_class):
    """ Check if class in list of REACTION CLASS
    """
    return rxn_class in _values(ReactionClass)


def reverse_reaction_class(rxn_class):
    """ determine the reverse of a reaction class
    """
    return REVERSE_REACTION_DCT.get(rxn_class, None)


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
