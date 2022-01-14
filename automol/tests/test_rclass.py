""" Test functionality that handles various reaction class
    desginations and descriptions
"""

import automol.par
from automol.par import ReactionClass


RCLASS1 = (ReactionClass.Typ.HYDROGEN_ABSTRACTION,
           ReactionClass.HIGHSPIN,
           False,
           False)
RCLASS2 = (ReactionClass.Typ.HYDROGEN_ABSTRACTION,
           ReactionClass.LOWSPIN,
           True,
           False)
RCLASS3 = (ReactionClass.Typ.HYDROGEN_ABSTRACTION,
           ReactionClass.HIGHSPIN,
           True,
           False)
RCLASS4 = (ReactionClass.Typ.ADDITION,
           ReactionClass.HIGHSPIN,
           False,
           False)
RCLASS5 = (ReactionClass.Typ.ADDITION,
           ReactionClass.LOWSPIN,
           True,
           False)
RCLASS6 = (ReactionClass.Typ.SUBSTITUTION,
           ReactionClass.LOWSPIN,
           True,
           False)
RCLASS7 = (ReactionClass.Typ.HYDROGEN_MIGRATION,
           '',
           False,
           False)
RCLASS8 = (ReactionClass.Typ.ELIMINATION,
           '',
           False,
           False)
RCLASS9 = (ReactionClass.Typ.ADDITION,
           None,
           False,
           True)


def test__build():
    """ test autmol.par.reaction_class_from_data
    """

    assert RCLASS1 == automol.par.reaction_class_from_data(
        ReactionClass.Typ.HYDROGEN_ABSTRACTION,
        ReactionClass.HIGHSPIN,
        False,
        False)
    assert RCLASS5 == automol.par.reaction_class_from_data(
        ReactionClass.Typ.ADDITION,
        ReactionClass.LOWSPIN,
        True,
        False)


def test__get():
    """ test automol.par.typ
        test automol.par.spin
        test automol.par.radrad
    """

    assert automol.par.typ(RCLASS1) == (
            ReactionClass.Typ.HYDROGEN_ABSTRACTION)
    assert automol.par.spin(RCLASS1) == (
            ReactionClass.HIGHSPIN)
    assert not automol.par.radrad(RCLASS1)
    assert not automol.par.isc(RCLASS1)

    assert automol.par.typ(RCLASS6) == (
            ReactionClass.Typ.SUBSTITUTION)
    assert automol.par.spin(RCLASS6) == (
            ReactionClass.LOWSPIN)
    assert automol.par.radrad(RCLASS6)
    assert not automol.par.isc(RCLASS1)

    assert automol.par.typ(RCLASS9) == (
            ReactionClass.Typ.ADDITION)
    assert automol.par.spin(RCLASS9) is None
    assert not automol.par.radrad(RCLASS9)
    assert automol.par.isc(RCLASS9)


def test__prop():
    """ test automol.par.is_high_spin
        test automol.par.is_low_spin
        test automol.par.has_nobarrier
    """

    assert automol.par.is_high_spin(RCLASS4)
    assert not automol.par.is_high_spin(RCLASS6)

    assert not automol.par.is_low_spin(RCLASS4)
    assert automol.par.is_low_spin(RCLASS6)

    assert not automol.par.is_radrad(RCLASS4)
    assert automol.par.is_radrad(RCLASS5)

    assert not automol.par.has_nobarrier(RCLASS1)
    assert automol.par.has_nobarrier(RCLASS2)
    assert not automol.par.has_nobarrier(RCLASS3)

    assert not automol.par.is_isc(RCLASS1)
    assert automol.par.is_isc(RCLASS9)


def test__need():
    """ test automol.par.need_spin_designation
        test automol.par.need_wells
    """

    assert automol.par.need_spin_designation(RCLASS1)
    assert automol.par.need_spin_designation(RCLASS4)
    assert not automol.par.need_spin_designation(RCLASS7)

    assert automol.par.need_wells(RCLASS1)
    assert automol.par.need_wells(RCLASS6)
    assert not automol.par.need_wells(RCLASS7)


def test__string():
    """ test automol.par.string
    """

    assert automol.par.string(RCLASS1) == (
        'high-spin hydrogen abstraction')
    assert automol.par.string(RCLASS2) == (
         'radical-radical low-spin hydrogen abstraction')
    assert automol.par.string(RCLASS9) == (
         'intersystem-crossing addition')
