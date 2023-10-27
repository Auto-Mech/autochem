""" Handle Full Reaction-Class designations for reactions that
    describe all meaningful attributes of the reaction required
    for electronic structure and kinetic calculations.
"""
import dataclasses
import enum


class ReactionClass(str, enum.Enum):
    """Reaction class names"""

    TRIVIAL = "trivial"
    # Unimolecular reactions
    HYDROGEN_MIGRATION = "hydrogen migration"
    BETA_SCISSION = "beta scission"
    HOMOLYT_SCISSION = "homolytic scission"
    RING_FORM_SCISSION = "ring forming scission"
    ELIMINATION = "elimination"
    # Bimolecular reactions
    HYDROGEN_ABSTRACTION = "hydrogen abstraction"
    ADDITION = "addition"
    INSERTION = "insertion"
    DOUBLE_INSERTION = "double_insertion"
    SUBSTITUTION = "substitution"

    def __str__(self):
        return self.value

    def __repr__(self):
        return repr(self.value)

    @classmethod
    def reverse(cls, value: str):
        """Get the class for the reverse of a reaction

        :param value: A reaction class
        :type value: str
        :return: The class of the reverse reaction
        :rtype: ReactionClass
        """
        reverse_dct = {
            cls.TRIVIAL: cls.TRIVIAL,
            # Unimolecular reactions
            cls.HYDROGEN_MIGRATION: cls.HYDROGEN_MIGRATION,
            cls.BETA_SCISSION: cls.ADDITION,
            cls.HOMOLYT_SCISSION: None,
            cls.RING_FORM_SCISSION: None,
            cls.ELIMINATION: cls.INSERTION,
            # Bimolecular reactions
            cls.HYDROGEN_ABSTRACTION: cls.HYDROGEN_ABSTRACTION,
            cls.ADDITION: cls.BETA_SCISSION,
            cls.INSERTION: cls.ELIMINATION,
            cls.DOUBLE_INSERTION: None,
            cls.SUBSTITUTION: cls.SUBSTITUTION,
        }
        return reverse_dct[cls(value)]

    @classmethod
    def is_reversible(cls, value: str) -> bool:
        """Is this reaction class reversible?

        :param value: A reaction class
        :type value: str
        :return: `True` if it is, `False` if it isn't
        :rtype: bool
        """
        return cls.reverse(value) is not None

    @classmethod
    def is_bimolecular(cls, value: str) -> bool:
        """Is this reaction class bimolecular?

        :param value: The reaction class
        :type value: str
        :return: `True` if it is, `False` if it isn't
        :rtype: bool
        """
        bimol_classes = (
            cls.HYDROGEN_ABSTRACTION,
            cls.ADDITION,
            cls.INSERTION,
            cls.DOUBLE_INSERTION,
            cls.SUBSTITUTION,
        )
        return cls(value) in bimol_classes

    @classmethod
    def requires_spin_designation(cls, value: str) -> bool:
        """Is this a reaction class that requires a spin designation?

        :param value: The reaction class
        :type value: str
        :return: `True` if it is, `False` if it isn't
        :rtype: bool
        """
        need_spin_classes = (
            cls.HYDROGEN_ABSTRACTION,  # AVC: Why is this in here??
            cls.ADDITION,
        )
        return cls(value) in need_spin_classes

    @classmethod
    def is_defined(cls, value: str) -> bool:
        """Is this reaction class defined?

        :param value: The reaction class
        :type value: str
        :return: `True` if it is, `False` if it isn't
        :rtype: bool
        """
        return str(value) in list(cls)


class ReactionSpin(str, enum.Enum):
    """reaction spin types"""

    LOW = "low-spin"
    HIGH = "high-spin"
    NONE = "unspecified spin"

    def __str__(self):
        return self.value

    def __repr__(self):
        return repr(self.value)

    @classmethod
    def is_defined(cls, value: str) -> bool:
        """Check whether a reaction spin type is defined

        :param value: The value to check for
        :type value: str
        :return: `True` if it does, `False` if it doesn't
        :rtype: bool
        """
        return value in list(cls)


@dataclasses.dataclass
class ReactionInfo:
    """General information about a reaction

    :param class_: The class name of the reaction
    :type class_: ReactionClass
    :param spin_: The spin-type of the reaction (low or high)
    :type spin_: ReactionSpin
    :param is_rad_rad: Whether this is a radical-radical reaction
    :type is_rad_rad: bool
    :param is_isc: Whether this is an intersystem crossing
    :type is_isc: bool
    """

    class_: ReactionClass
    spin: ReactionSpin = ReactionSpin.NONE
    is_rad_rad: bool = False
    is_isc: bool = False

    def __str__(self) -> str:
        """Generate a string representation of the reaction information"""
        parts = []

        if self.is_radical_radical():
            parts.append("radical-radical")

        if self.is_intersystem_crossing():
            parts.append("intersystem-crossing")

        if not self.has_no_spin_designation():
            parts.append(self.reaction_spin())

        parts.append(self.reaction_class())
        return " ".join(map(str, parts))

    def __repr__(self) -> str:
        """Generate a string representation of the reaction information"""
        return str(self)

    def string(self) -> str:
        """Get a string representation of the reaction"""
        return str(self)

    def reaction_class(self) -> ReactionClass:
        """Get the reaction class name"""
        return ReactionClass(self.class_)

    def reaction_spin(self) -> ReactionSpin:
        """Get the reaction spin type"""
        return ReactionSpin(self.spin)

    def is_radical_radical(self) -> bool:
        """Is this a radical radical reaction?"""
        return self.is_rad_rad

    def is_intersystem_crossing(self) -> bool:
        """Is this an intersystem crossing?"""
        return self.is_isc

    def is_barrierless(self) -> bool:
        """Is this a barrierless reaction?"""
        return self.is_radical_radical() and not self.is_high_spin()

    def is_low_spin(self) -> bool:
        """Is this a low-spin reaction?"""
        return self.reaction_spin() == ReactionSpin.LOW

    def is_high_spin(self) -> bool:
        """Is this a high-spin reaction?"""
        return self.reaction_spin() == ReactionSpin.HIGH

    def has_no_spin_designation(self) -> bool:
        """Is this the only possible spin-state?"""
        return self.reaction_spin() == ReactionSpin.NONE

    def requires_spin_designation(self) -> bool:
        """Is a spin designation required for this reaction?"""
        spin_req_classes = (
            ReactionClass.HYDROGEN_ABSTRACTION,  # AVC: Why is this in here??
            ReactionClass.ADDITION,
        )
        return (
            self.reaction_class() in spin_req_classes
            and not self.is_intersystem_crossing()
        )

    def requires_well_description(self) -> bool:
        """Determine if a reaction is appropriately described by the presence of
        entrance- or exit-channel van der Waals wells
        """
        well_classes = (
            ReactionClass.HYDROGEN_ABSTRACTION,
            ReactionClass.SUBSTITUTION,
        )
        return self.reaction_class() in well_classes
