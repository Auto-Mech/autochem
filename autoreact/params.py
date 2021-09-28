""" Handle RxnParams objects
"""

import numpy

class RxnParams:
    """ Store and manipulate parameters that correspond to
        various functional form expressions used to represent
        reaction rate constants for some range of temperature and pressure
    """

    def __init__(self, arr_dct=None, plog_dct=None, cheb_dct=None,
                 troe_dct=None, lind_dct=None):

        self.arr = None
        self.plog = None
        self.cheb = None
        self.troe = None
        self.lind = None

        if arr_dct is not None:
            self.set_arr(arr_dct)
        if plog_dct is not None:
            self.set_plog(plog_dct)
        if cheb_dct is not None:
            self.set_cheb(cheb_dct)
        if troe_dct is not None:
            pass

    # Parameter setting functions
    def set_arr(self, arr_dct):
        """ Sets Arrhenius parameters

        """

        for key in arr_dct.keys():
            assert key in ('arr_tuples', 'arr_collid'), (
                f'arr_dct keys should be "arr_tuples" or "arr_collid", not "{key}"')        

        arr_tuples = arr_dct.get('arr_tuples')
        self.arr_collid = arr_dct.get('arr_collid')
        self.check_arr(arr_tuples)
        if self.arr is None:
            self.arr = tuple(tuple(arr_tuple) for arr_tuple in arr_tuples)
        else:
            self.arr += arr_tuples


    def set_plog(self, plog_dct):
        """ Sets PLOG parameters.

            :param plog_dct: Arrhenius fitting parameters at pressures
            :type plog_dct: dict[float: tuple(tuple(float))]
        """

        # Check for the bad_dup flag; if found, remove from plog_dct and flag
        self.bad_dup = False
        if 'bad_dup' in plog_dct:
            plog_dct.pop('bad_dup')
            self.bad_dup = True

        # Add the PLOG parameters to the RxnParams object
        self.check_plog(plog_dct)
        if self.plog is None:
            self.plog = {pressure: arr_tuples
                         for pressure, arr_tuples in plog_dct.items()
                         if pressure != 'high'}
        else:
            print("Only allowed one PLOG set")


    def set_cheb(self, cheb_dct):
        """ Sets Chebyshev parameters. 
        """

        self.check_cheb(cheb_dct)

        # Fill the one_atm_arr field if it's None (this may be redundant)
        if cheb_dct.get('one_atm_arr') is None:
            cheb_dct['one_atm_arr'] = None

        if self.cheb is None:
            self.cheb = cheb_dct
        else:
            print("Only allowed one Chebyshev set")


    def set_troe(self, troe_dct): 
        """ Sets Troe parameters. """

        RxnParams.check_troe_lind(highp_arr, lowp_arr,
                                  troe_params=troe_params, alpha=alpha)
        if self.troe is None:
            self.troe = troe_dct
        else:
            print("Only allowed one Troe set")


    def set_lind(self, highp_arr, lowp_arr, colliders=None):
        """ Sets Lindemann parameters. """

        RxnParams.check_troe_lind(highp_arr, lowp_arr)
        if self.lind is None:
            self.lind = {
                'highp_arr': highp_arr,
                'lowp_arr': lowp_arr,
            }
            self.colliders = colliders
        else:
            print("Only allowed one Lindemann set")
            

    def check_arr(self, arr_tuples):
        """ Ensures that input Arrhenius parameters have correct form:
            ((a1, n1, e1), (a2, n2, e2), ...)

            :param arr_tuples: Arrhenius parameters
            :type arr_tuples: tuple(tuple(float))
        """
        assert isinstance(arr_tuples, (list, tuple)), (
            'arr_tuples should be a list-of-lists or tuple-of-tuples')
        for arr_tuple in arr_tuples:
            assert isinstance(arr_tuple, (list, tuple)), (
                'arr_tuple should be a list or tuple')
            assert len(arr_tuple) == 3, (
                'length of each arr_tuple should be three')
            assert all(isinstance(x, (int, float)) for x in arr_tuple), (
                'each parameter in the param set should be a float')


    def check_plog(self, plog_dct):
        """ Ensures that input PLOG dictionary has the correct form:
            {'high'/pressure: ((a1, n1, e1), (a2, n2, e2), ...)}

            :param plog_dct: Arrhenius fitting parameters at pressures
            :type plog_dct: dict[float: tuple(tuple(float))]
        """
        assert isinstance(plog_dct, dict)
        for pressure, arr_tuples in plog_dct.items():
            assert pressure == 'high' or isinstance(pressure, float)
            self.check_arr(arr_tuples)


    def check_cheb(self, cheb_dct):
        """ Ensure that the Chebyshev parameters are of the proper form
        """

        tlim = cheb_dct.get('tlim')
        plim = cheb_dct.get('plim')
        alpha= cheb_dct.get('alpha')
        one_atm_arr = cheb_dct.get('one_atm_arr')

        if one_atm_arr is not None:  # the one atm rates are not required
            self.check_arr(one_atm_arr)        
        assert isinstance(tlim, (list, tuple))
        assert isinstance(plim, (list, tuple))
        assert all(isinstance(x, (int, float)) for x in tlim)
        assert all(isinstance(x, (int, float)) for x in plim)
        assert isinstance(alpha, (list, tuple, numpy.ndarray))
        for row in alpha:
            assert isinstance(row, (list, tuple, numpy.ndarray))
            assert all(isinstance(x, (int, float)) for x in row)


    def check_troe_lind(self, param_dct):
        """ Ensures that the Troe or Lindemann parameters are of proper form
        """

        highp_arr = param_dct['highp_arr']
        lowp_arr = param_dct['lowp_arr']
        troe_params = param_dct.get('troe_params')

        check_arr(highp_arr)
        check_arr(lowp_arr)
        if troe_params is not None:
            assert len(troe_params) in (3, 4)
            assert all(isinstance(x, float) for x in troe_params)


    def existing_parameters(self):
        """ Return a list of all the types of parameters contained
            in the object
        """

        _existing = ()
        if self.arr is not None:
            _existing += ('arr',)
        if self.plog is not None:
            _existing += ('plog',)
        if self.cheb is not None:
            _existing += ('cheb',)
        if self.troe is not None:
            _existing += ('troe',)
        if self.lind is not None:
            _existing += ('lind',)

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
            elif 'cheb' in _existing_params:
                _write = 'cheb'
            elif 'troe' in _existing_params:
                _write = 'troe'
            elif 'lind' in _existing_params:
                _write = 'lind'
            elif 'arr' in _existing_params:
                _write = 'arr'
        else:
            _write = None

        return _write
