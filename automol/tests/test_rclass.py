""" Test functionality that handles various reaction class
    desginations and descriptions
"""
import automol
from automol import ReactionClass, ReactionInfo, ReactionSpin

RINFO1 = ReactionInfo(
    class_=ReactionClass.HYDROGEN_ABSTRACTION,
    spin=ReactionSpin.HIGH,
    is_rad_rad=False,
    is_isc=False,
)
RINFO2 = ReactionInfo(
    class_=ReactionClass.HYDROGEN_ABSTRACTION,
    spin=ReactionSpin.LOW,
    is_rad_rad=True,
    is_isc=False,
)
RINFO3 = ReactionInfo(
    class_=ReactionClass.HYDROGEN_ABSTRACTION,
    spin=ReactionSpin.HIGH,
    is_rad_rad=True,
    is_isc=False,
)
RINFO4 = ReactionInfo(
    class_=ReactionClass.ADDITION,
    spin=ReactionSpin.HIGH,
    is_rad_rad=False,
    is_isc=False,
)
RINFO5 = ReactionInfo(
    class_=ReactionClass.ADDITION,
    spin=ReactionSpin.LOW,
    is_rad_rad=True,
    is_isc=False,
)
RINFO6 = ReactionInfo(
    class_=ReactionClass.SUBSTITUTION,
    spin=ReactionSpin.LOW,
    is_rad_rad=True,
    is_isc=False,
)
RINFO7 = ReactionInfo(
    class_=ReactionClass.HYDROGEN_MIGRATION,
    is_rad_rad=False,
    is_isc=False,
)
RINFO8 = ReactionInfo(
    class_=ReactionClass.ELIMINATION,
    is_rad_rad=False,
    is_isc=False,
)
RINFO9 = ReactionInfo(
    class_=ReactionClass.ADDITION,
    is_rad_rad=False,
    is_isc=True,
)


def test__build():
    """test autmol.par.reaction_class_from_data"""

    assert RINFO1 == ReactionInfo(
        class_=ReactionClass.HYDROGEN_ABSTRACTION,
        spin=ReactionSpin.HIGH,
        is_rad_rad=False,
        is_isc=False,
    )
    assert RINFO5 == ReactionInfo(
        class_=ReactionClass.ADDITION,
        spin=ReactionSpin.LOW,
        is_rad_rad=True,
        is_isc=False,
    )


def test__get():
    """test automol.par.typ
    test automol.par.spin
    test automol.par.radrad
    """

    assert ReactionInfo.reaction_class(RINFO1) == ReactionClass.HYDROGEN_ABSTRACTION
    assert ReactionInfo.reaction_spin(RINFO1) == ReactionSpin.HIGH
    assert not ReactionInfo.is_radical_radical(RINFO1)
    assert not ReactionInfo.is_intersystem_crossing(RINFO1)

    assert ReactionInfo.reaction_class(RINFO6) == ReactionClass.SUBSTITUTION
    assert ReactionInfo.reaction_spin(RINFO6) == ReactionSpin.LOW
    assert ReactionInfo.is_radical_radical(RINFO6)
    assert not ReactionInfo.is_intersystem_crossing(RINFO1)

    assert ReactionInfo.reaction_class(RINFO9) == ReactionClass.ADDITION
    assert ReactionInfo.reaction_spin(RINFO9) is ReactionSpin.NONE
    assert not ReactionInfo.is_radical_radical(RINFO9)
    assert ReactionInfo.is_intersystem_crossing(RINFO9)


def test__prop():
    """test automol.par.is_high_spin
    test automol.par.is_low_spin
    test automol.par.has_nobarrier
    """

    assert ReactionInfo.is_high_spin(RINFO4)
    assert not ReactionInfo.is_high_spin(RINFO6)

    assert not ReactionInfo.is_low_spin(RINFO4)
    assert ReactionInfo.is_low_spin(RINFO6)

    assert not ReactionInfo.is_radical_radical(RINFO4)
    assert ReactionInfo.is_radical_radical(RINFO5)

    assert not ReactionInfo.is_barrierless(RINFO1)
    assert ReactionInfo.is_barrierless(RINFO2)
    assert not ReactionInfo.is_barrierless(RINFO3)

    assert not ReactionInfo.is_intersystem_crossing(RINFO1)
    assert ReactionInfo.is_intersystem_crossing(RINFO9)


def test__need():
    """test automol.par.need_spin_designation
    test automol.par.need_wells
    """

    assert automol.ReactionInfo.requires_spin_designation(RINFO1)
    assert automol.ReactionInfo.requires_spin_designation(RINFO4)
    assert not automol.ReactionInfo.requires_spin_designation(RINFO7)

    assert automol.ReactionInfo.requires_well_description(RINFO1)
    assert automol.ReactionInfo.requires_well_description(RINFO6)
    assert not automol.ReactionInfo.requires_well_description(RINFO7)


def test__string():
    """test automol.par.string"""

    assert ReactionInfo.string(RINFO1) == ("high-spin hydrogen abstraction")
    assert ReactionInfo.string(RINFO2) == (
        "radical-radical low-spin hydrogen abstraction"
    )
    assert ReactionInfo.string(RINFO9) == ("intersystem-crossing addition")


if __name__ == '__main__':
    test__prop()
