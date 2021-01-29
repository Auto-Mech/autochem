

def low_spin_multiplicity(ich):
    """ Guess spin multiplicity based on the number of electrons.

        :param ich: InChI string
        :type ich: str
        :rtype: int
    """

    fml = automol.inchi.formula(ich)
    nelec = automol.formula.electron_count(fml)

    if (nelec % 2) == 0:
        mult = 1
    else:
        mult = 2

    return mult
