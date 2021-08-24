""" Handle RxnParams objects
"""


class RxnParams:
    """ Store and manipulate parameters that correspond to
        various functional form expressions used to represent
        reaction rate constants for some range of temperature and pressure
    """

    def __init__(self):
        self.arrhenius = None
        self.plog = None
        self.chebyshev = None
        self.troe = None
        self.lindemann = None
        self.colliders = None

    # Parameter setting functions
    def set_arr(self, arr_params):
        """ Sets Arrhenius parameters

            :param arr_params: Arrhenius parameters
            :type arr_params: list [A, n, Ea]
        """

        RxnParams.check_arrhenius(arr_params)
        if self.arrhenius is None:
            self.arrhenius = tuple(tuple(par) for par in arr_params)
        else:
            self.arrhenius += arr_params

    def set_plog(self, plog_dct):
        """ Sets PLog parameters.

            :param plog_dct: Arrhenius fitting parameters at pressures
            :type plog_dct: dict[float: tuple(tuple(float))]
        """

        RxnParams.check_plog(plog_dct)
        if self.plog is None:
            self.plog = {pressure: params
                         for pressure, params in plog_dct.items()
                         if pressure != 'high'}
        else:
            print("Only allowed one PLog set")
            # for pressure in plog_dct.keys():
            #     if pressure == 'high':
            #         self_pressure = 'high'
            #     else:
            #         pass
            #         # for
            #     if pressure in self.plog:
            #         self.plog[self_pressure] += plog_dct[pressure]
            #     else:
            #         self.plog[self_pressure] = plog_dct[pressure]

    def set_chebyshev(self, highp_arr, lowp_arr, tlim, plim, alpha):
        """ Sets Chebyshev parameters. """

        RxnParams.check_chebyshev(highp_arr, lowp_arr, tlim, plim, alpha)
        if self.chebyshev is None:
            self.chebyshev = {
                'highp_arr': highp_arr,
                'lowp_arr': lowp_arr,
                'tlim': tlim,
                'plim': plim,
                'alpha': alpha
            }
        else:
            print("Only allowed one Chebyshev set")

    def set_troe(self, highp_arr, lowp_arr,
                 alpha, troe_params, colliders=None):
        """ Sets Troe parameters. """

        RxnParams.check_troe_lind(highp_arr, lowp_arr,
                                  troe_params=troe_params, alpha=alpha)
        if self.troe is None:
            self.troe = {
                'highp_arr': highp_arr,
                'lowp_arr': lowp_arr,
                'alpha': alpha,
                'troe_params': troe_params
            }
            self.colliders = colliders
        else:
            print("Only allowed one Troe set")

    def set_lindemann(self, highp_arr, lowp_arr, colliders=None):
        """ Sets Lindemann parameters. """

        RxnParams.check_troe_lind(highp_arr, lowp_arr)
        if self.lindemann is None:
            self.lindemann = {
                'highp_arr': highp_arr,
                'lowp_arr': lowp_arr,
            }
            self.colliders = colliders
        else:
            print("Only allowed one Lindemann set")

    @staticmethod
    def check_arrhenius(arr_params):
        """ Ensure that input Arrhenius parameters have correct form:
            ((a1, n1, e1), (a2, n2, e2), ...)

            :param arr_params: arrhenius parameters
            :type arr_params: tuple(tuple(float))
        """
        assert isinstance(arr_params, (list, tuple)), (
            'arr_params should be a list-of-lists or tuple-of-tuples')
        for params in arr_params:
            assert isinstance(params, (list, tuple)), (
                'params should be a list-of-lists or tuple-of-tuples')
            assert len(params) == 3, (
                'length of each params should be three')
            assert all(isinstance(x, (int, float)) for x in params), (
                'each parameter in the param set should be a float')

    @staticmethod
    def check_plog(plog_dct):
        """ Ensure that input PLog dictionary has the correct form:
            {'high'/pressure: ((a1, n1, e1), (a2, n2, e2), ...)}

            :param plog_dct: Arrhenius fitting parameters at pressures
            :type plog_dct: dict[float: tuple(tuple(float))]
        """
        assert isinstance(plog_dct, dict)
        for key, val in plog_dct.items():
            assert key == 'high' or isinstance(key, float)
            RxnParams.check_arrhenius(val)

    @staticmethod
    def check_chebyshev(highp_arr, lowp_arr, tlim, plim, alpha):
        """ Ensure that the Chebyshev parameters are of the proper form
        """
        RxnParams.check_arrhenius(highp_arr)
        RxnParams.check_arrhenius(lowp_arr)
        assert isinstance(tlim, (list, tuple))
        assert isinstance(plim, (list, tuple))
        assert all(isinstance(x, (int, float)) for x in tlim)
        assert all(isinstance(x, (int, float)) for x in plim)
        assert isinstance(alpha, (list, tuple))
        for row in alpha:
            assert isinstance(row, (list, tuple))
            assert all(isinstance(x, (int, float)) for x in row)

    @staticmethod
    def check_troe_lind(highp_arr, lowp_arr, alpha=None, troe_params=None):
        """ Enusre that the Troe or Lindemann parameters are of proper form
        """
        RxnParams.check_arrhenius(highp_arr)
        RxnParams.check_arrhenius(lowp_arr)
        if alpha is not None:
            isinstance(alpha, float)
        if troe_params is not None:
            assert len(troe_params) in (3, 4)
            assert all(isinstance(x, float) for x in troe_params)

    def existing_parameters(self):
        """ Return a list of all the types of parameters contained
            in the object
        """

        _existing = ()
        if self.arrhenius is not None:
            _existing += ('arrhenius',)
        if self.plog is not None:
            _existing += ('plog',)
        if self.chebyshev is not None:
            _existing += ('chebyshev',)
        if self.troe is not None:
            _existing += ('troe',)
        if self.lindemann is not None:
            _existing += ('lindemann',)

        return _existing

    # Function to choose the best functional expression
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
                _write = 'lindemann'
            elif 'arrhenius' in _existing_params:
                _write = 'arrhenius'
        else:
            _write = None

        return _write
