"""
  Build class definition
"""

class ZMAT():
    """ Z-Matrix class
    """

    def __init__(self, zma):
        self.zma = zma
        self.frm_bnd_keys = frozenset({})
        self.frm_brk_keys = frozenset({})
        self.tors_name_keys = frozenset({})
        self.const_name_keys = frozenset({})
