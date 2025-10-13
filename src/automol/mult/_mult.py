"""Calculate spins (2S) from spin multiplicities."""


def spin(mult: int) -> int:
    """Calculate the spin (2Ms) from using the spin multiplicity
    which is equivalent to the number of unpaired electrons.

    :param mult: Multiplicity
    :return: Number of unpaired electrons
    """
    return mult - 1
