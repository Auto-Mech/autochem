""" Handle Full Reaction-Class designations for reactions that
    describe all meaningful attributes of the reaction required
    for electronic structure and kinetic calculations.

    rclass = (reaction type, high/low spin, radical-radical)
             (str, str, bool)

    reaction type: corresponds to labels set in ReactionClass.Typ

"""

import inspect


class ReactionClass:
    """ Full Description of a Reaction Class
    """
    class Typ:
        """ Names of supported reaction classes
        """
        TRIVIAL = 'trivial'
        # Unimolecular reactions
        HYDROGEN_MIGRATION = 'hydrogen migration'
        BETA_SCISSION = 'beta scission'
        HOMOLYT_SCISSION = 'homolytic scission'
        RING_FORM_SCISSION = 'ring forming scission'
        ELIMINATION = 'elimination'
        # Bimolecular reactions
        HYDROGEN_ABSTRACTION = 'hydrogen abstraction'
        ADDITION = 'addition'
        INSERTION = 'insertion'
        SUBSTITUTION = 'substitution'

    class MechTyp:
        """ Names of roles reaction serves in a mechanism
        """
        PROPAGATION = 'propagation'
        TERMINATION = 'termination'
        BRANCHING = 'branching'
        LUMPED = 'lumped'

    # class Prop:
    RADRAD = 'radical-radical'
    LOWSPIN = 'low-spin'
    HIGHSPIN = 'high-spin'
    ISC = 'intersystem-crossing'


REVERSE_REACTION_DCT = {
    ReactionClass.Typ.HYDROGEN_MIGRATION:
        ReactionClass.Typ.HYDROGEN_MIGRATION,
    ReactionClass.Typ.HYDROGEN_ABSTRACTION:
        ReactionClass.Typ.HYDROGEN_ABSTRACTION,
    ReactionClass.Typ.ADDITION: ReactionClass.Typ.BETA_SCISSION,
    ReactionClass.Typ.BETA_SCISSION: ReactionClass.Typ.ADDITION,
    ReactionClass.Typ.ELIMINATION: ReactionClass.Typ.INSERTION,
    ReactionClass.Typ.INSERTION: ReactionClass.Typ.ELIMINATION,
    ReactionClass.Typ.SUBSTITUTION: None
}

BIMOL_REACTIONS = [
    ReactionClass.Typ.HYDROGEN_ABSTRACTION,
    ReactionClass.Typ.ADDITION,
    ReactionClass.Typ.INSERTION,
    ReactionClass.Typ.SUBSTITUTION
]


# Builders
def reaction_class_from_data(class_typ, class_spin,
                             class_radrad, class_isc):
    """ Build a full-class description including the following useful
        descriptors of the reaction class:

            typ: type of reaction (e.g., abstraction, addition, migration)
            spin: whether the reaction occurs on a low- or high-spin state
            radrad: whether the reactants or products are two-or-more radicals

        :param class_typ: reaction type designation
        :type class_typ: str
        :param class_spin: spin type designation
        :type class_spin: str
        :param class_radrad: radical-radical type designation
        :type class_radrad: bool
        :param class_isc: intersystem-crossing designation
        :type class_isc: bool
        :rtype: (str, str, bool, bool)
    """
    return (class_typ, class_spin, class_radrad, class_isc)


def string(rxn_class):
    """ Write class to string

        :param rxn_class: reaction class including type, spin, radrad
        :type rxn_class:
        :rtype: str
    """
    spin_str = '' if spin(rxn_class) is None else spin(rxn_class)
    radrad_str = '' if not radrad(rxn_class) else 'radical-radical'
    isc_str = '' if not isc(rxn_class) else 'intersystem-crossing'
    cls_str = typ(rxn_class)

    out_str = ''
    for str_ in (isc_str, radrad_str, spin_str, cls_str):
        if str_:
            out_str += '{} '.format(str_)

    return out_str.strip()


# Checks for building/using class
def need_spin_designation(rxn_class):
    """ Determine if a spin-state string designation in the full reaction
        class description based on the class typ

        :rtype: bool
    """
    need_spins = (
        ReactionClass.Typ.ADDITION,
        ReactionClass.Typ.HYDROGEN_ABSTRACTION
    )
    return (typ(rxn_class) in need_spins) and not isc(rxn_class)


def need_wells(rxn_class):
    """ Determine if a reaction is appropriately described
        by the presence of entrance- or exit-channel van der Waals
        wells.
    """
    _need_wells = (
        ReactionClass.Typ.HYDROGEN_ABSTRACTION,
        ReactionClass.Typ.SUBSTITUTION
    )
    return typ(rxn_class) in _need_wells


# Get designations from the reaction class
def typ(rxn_class):
    """ Get the reaction type designation from the reaction class.

        :param rxn_class: reaction class including type, spin, radrad
        :type rxn_class:
        :rtype: bool
    """
    return rxn_class[0]


def spin(rxn_class):
    """ Get the spin designation from the reaction class.

        :param rxn_class: reaction class including type, spin, radrad
        :type rxn_class:
        :rtype: bool
    """
    return rxn_class[1]


def radrad(rxn_class):
    """ Get the radical-radical designation from the reaction class.

        :param rxn_class: reaction class including type, spin, radrad
        :type rxn_class:
        :rtype: bool
    """
    return rxn_class[2]


def isc(rxn_class):
    """ Get the intersystem crossing designation from the reaction class.

        :param rxn_class: reaction class including type, spin, radrad
        :type rxn_class:
        :rtype: bool
    """
    return rxn_class[3]


def is_high_spin(rxn_class):
    """ Return boolean for whether a reaction is a high-spin reaction
        based on the spin designation in its class description.

        :param rxn_class: reaction class including type, spin, radrad
        :type rxn_class:
        :rtype: bool
    """
    return bool(ReactionClass.HIGHSPIN in spin(rxn_class))


def is_low_spin(rxn_class):
    """ Return boolean for whether a reaction is a low-spin reaction
        based on the spin designation in its class description.

        :param rxn_class: reaction class including type, spin, radrad
        :type rxn_class:
        :rtype: bool
    """
    return bool(ReactionClass.LOWSPIN in spin(rxn_class))


def is_radrad(rxn_class):
    """ Return boolean for whether a reaction is a radical-radical reaction
        based on the spin designation in its class description.

        :param rxn_class: reaction class including type, spin, radrad
        :type rxn_class:
        :rtype: bool
    """
    return radrad(rxn_class)


def has_nobarrier(rxn_class):
    """ Determine if reaction is (likely) to be barrierless.

        :param rxn_class: reaction class including type, spin, radrad
        :type rxn_class:
        :rtype: bool
    """
    return is_radrad(rxn_class) and not is_high_spin(rxn_class)


def is_isc(rxn_class):
    """ Return boolean for whether a reaction is an
        intrsystem crossing reaction

        :param rxn_class: reaction class including type, spin, radrad
        :type rxn_class:
        :rtype: bool
    """
    return isc(rxn_class)


# Handle Reaction Class IDs
def is_reaction_class(rxn_class):
    """ Check if class in list of REACTION CLASS
    """
    return rxn_class in _values(ReactionClass.Typ)


def reverse_reaction_class(rxn_class):
    """ determine the reverse of a reaction class
    """
    return REVERSE_REACTION_DCT.get(rxn_class, None)


# Functions for just the reactant type
def isbimol(rxn_typ):
    """ Determine if a type is a bimolecular reaction
    """
    return rxn_typ in BIMOL_REACTIONS


# Other Utility Functions
def _values(cls):
    """ list the values of a parameter class
    """
    assert inspect.isclass(cls)
    vals = tuple(val for val in _public_attributes(cls)
                 if not inspect.isclass(val))
    return vals


def _public_attributes(cls):
    return tuple(val for name, val in
                 inspect.getmembers(cls, lambda x: not inspect.isroutine(x))
                 if not name.startswith('_') and not inspect.isfunction(val))
