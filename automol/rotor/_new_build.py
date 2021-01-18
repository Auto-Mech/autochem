

class Rotor:
    """ Describes a Rotor
    """
    
    def __init__(self, names, zma, group_dct=None):
        """ constructor
        """
        self.names = names
        self.zma = zma

        torsions = ()
        for name in names:




class Torsion:
    """ Describes a torsion
    """

    def __init__(self, name, zma, group_dct=None):
        """ constructor
        """
        self.name = name
        self.zma = zma
        if group_dct is not None:
            group_dct = get(Parameter.ZMA, None)
            group_dct = get(Parameter.ZMA, None)
            group_dct = get(Parameter.ZMA, None)
            group_dct = get(Parameter.ZMA, None)
        else:
            self.symmetry = symmetry_number(self)
            self.span = span(self)
            self.axis, self.groups = rotational_groups(self)

    def symmetry_number(self):
        """
        """
        return 1.0


    def span(self):
        """
        """
        return (2.0 * numpy.pi) / symmetry_number(self)


    def rotational_groups(self):
        """
        """
        rot_group_dct = { (key1, key2): (group1, group2) }
        
