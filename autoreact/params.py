""" Handle RxnParams objects
"""


class RxnParams:
    """ Store and manipulate parameters that correspond to
        various functional form expressions used to represent
        reaction rate constants for some range of temperature and pressure
    """

    def __init__(self):
        self.arr = None
        self.plog = None
        self.cheb = None
        self.troe = None
        self.lind = None
        self.collid = None

    def set_arr(self, arr_params, arr_params2=None):
        """ Sets Arrhenius parameters

            :param arr_params: Arrhenius parameters
            :type arr_params: list [A, n, Ea]
            :param arr_params2: optional second set of Arrhenius parameters
            :type arr_params2: list [A, n, Ea]
        """

        assert isinstance(arr_params, list), 'arr_params should be a list'
        assert len(arr_params) == 3, (
            'Length of arr_params should be 3 but'
            ' is {len(arr_params)}')

        # Add the arr_params to the RxnParams object
        if self.arr is None:
            self.arr = [arr_params]
        else:
            self.arr.append(arr_params)

        # If a second set of Arrhenius parameters was input
        if arr_params2 is not None:
            assert isinstance(arr_params2, list), (
                'arr_params2 should be a list')
            assert len(arr_params2) == 3, (
                'Length of arr_params2 should be 3'
                ' but is {len(arr_params2)}')
            self.arr.append(arr_params2)

    def set_plog():
        """ Sets PLog parameters. """
        return NotImplementedError

    def set_cheb():
        """ Sets Chebyshev parameters. """
        return NotImplementedError

    def set_troe():
        """ Sets Troe parameters. """
        return NotImplementedError

    def existing_parameters(self):
        """ Return a list of all the types of parameters contained
            in the object
        """

        _existing = ()
        if self.arr is not None:
            _existing += ('arrhenius',)
        if self.plog is not None:
            _existing += ('plog',)
        if self.cheb is not None:
            _existing += ('chebyshev',)
        if self.troe is not None:
            _existing += ('troe',)
        if self.lind is not None:
            _existing += ('lindemann',)

        return _existing

    def form_to_write(self):
        """ Function that identifies which functional form should be
            written based on what parameters exist
        """

        _existing_params = self.existing_parameters()

        if _existing_params:
            if 'plog' in _existing_params:
                _write = 'plog'
            elif 'chebyshev' in _existing_params:
                _write = 'chebyshev'
            elif 'troe' in _existing_params:
                _write = 'troe'
            elif 'lindemann' in _existing_params:
                _write = 'linedmann'
            elif 'arrhenius' in _existing_params:
                _write = 'arrhenius'
        else:
            _write = None

        return _write
