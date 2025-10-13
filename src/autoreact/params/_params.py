""" Build and Manipulate RxnParams objects

TODO: After confirming that one reaction can't simultaneously have multiple
distinct parametriztaions, this class should be redesigned.

class RateParametrization():

    def __init__(self, type_, dict_):
        '''
            :param type_: The parametrization type. Options: 'arr', 'plog',
                'che', 'troe', 'lind'.
            :param dct: A dictionary of parameters corresponding to the
                parametrization type.
        '''
        self.type = type_
        self.dict = dict_
"""
import copy
import ast
import numpy
import pyparsing as pp


class RxnParams:
    """ Stores and manipulates parameters that correspond to
        various functional form expressions used to represent
        reaction rate constants for some range of temperature and pressure

        This is a horrendously designed class. Needs to be refactored.
    """

    def __init__(self, arr_dct=None, plog_dct=None, cheb_dct=None,
                 troe_dct=None, lind_dct=None):
        """ Creates an instance of RxnParams given an input dct

            :param arr_dct: dct describing Arrhenius parameters
            :type arr_dct: dct {'arr_tuples': (arr_tuple1, ...),
                'arr_collid': {spc1: eff1, spc2: ...}
            :param plog_dct: dct describing PLOG parameters
            :type plog_dct: dct {pressure1: arr_tuples1, pressure2: ...}
            :param cheb_dct: dct describing Chebyshev parameters
            :type cheb_dct: dct {'tlim': (tmin, tmax), 'plim': (pmin, pmax),
                'alpha': matrix of Cheb coefficients}
            :param troe_dct: dct describing Troe parameters
            :type troe_dct: dct {'highp_arr': highp_arr_tuples,
                'lowp_arr': lowp_arr_tuples,
                'troe_params': [alpha, T***, T*, T**],
                'collid': {spc1: eff1, spc2: ...}}
            :param lind_dct: dct describing Lindemann parameters
            :type lind_dct: dct {'highp_arr': highp_arr_tuples,
                'lowp_arr': lowp_arr_tuples, 'collid': {spc1: eff1, spc2: ...}}
        """

        self.arr = None
        self.plog = None
        self.cheb = None
        self.troe = None
        self.lind = None

        self.plog_dups = None  # no arr_dups since included in arr atribute
        self.cheb_dups = None
        self.troe_dups = None
        self.lind_dups = None

        if arr_dct is not None:
            self.set_arr(arr_dct)
        if plog_dct is not None:
            self.set_plog(plog_dct)
        if cheb_dct is not None:
            self.set_cheb(cheb_dct)
        if troe_dct is not None:
            self.set_troe(troe_dct)
        if lind_dct is not None:
            self.set_lind(lind_dct)

    @classmethod
    def from_string(cls, string):
        """ re-generated an instance from a string representation
        """
        word = pp.Word(pp.alphas)
        rest = pp.Regex('.*')
        parser = (word + '(' + word('type') + ',' + rest('params'))
        parse_dct = parser.parseString(string[:-1]).asDict()
        type_ = parse_dct['type']
        if type_ != 'arr':
            raise NotImplementedError(f"Parameter type '{type_}' not "
                                      f"yet implemented.")
        params = ast.literal_eval(parse_dct['params'])
        return cls(arr_dct={'arr_tuples': params})

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        """ Get a string representation of the parametrization
        """
        type_, dct = (('arr', (self.arr, self.arr_collid)) if self.arr is not None else
                      ('plog', self.plog) if self.plog is not None else
                      ('cheb', self.cheb) if self.cheb is not None else
                      ('troe', self.troe) if self.troe is not None else
                      ('lind', self.lind) if self.lind is not None else
                      (None, None))
        assert type_ is not None, 'No parametrization was assigned!'
        return f"RxnParams({type_}, {str(dct)})"

    def __radd__(self, other):
        ret = self
        if other != 0:
            ret = self.__add__(other)
        return ret

    def __add__(self, other):
        params = copy.deepcopy(self)
        params.combine_objects(other)
        return params

    def __scale__(self, scale_A, scale_n, scale_EA):
        return scale_factor(self, factor_A=scale_A, factor_n=scale_n, factor_EA=scale_EA)
    
    def __mul__(self, number):
        return scale_factor(self, factor_A=number)
    
    def __rmul__(self, number):
        return scale_factor(self, factor_A=number)

    def __truediv__(self, number):
        return scale_factor(self, factor_A=1/number)

    def __div__(self, number):
        return scale_factor(self, factor_A=1/number)

    def set_arr(self, arr_dct):
        """ Sets Arrhenius parameters

            :param arr_dct: dct describing Arrhenius parameters
            :type arr_dct: dct {'arr_tuples': (arr_tuple1, ...),
                'arr_collid': {spc1: eff1, spc2: ...}
        """

        # Run some checks
        for key in arr_dct.keys():
            assert key in ('arr_tuples', 'arr_collid'), (
                f'arr_dct keys should be "arr_tuples" or "arr_collid"'
                f', not "{key}"')
        arr_tuples = arr_dct.get('arr_tuples')
        collid = arr_dct.get('arr_collid')
        self.check_arr(arr_tuples)
        if collid is not None:
            self.check_collid(collid)

        # Add the arr params to the RxnParams object
        if self.arr is None:
            self.arr = tuple(tuple(arr_tuple) for arr_tuple in arr_tuples)
        else:
            self.arr += arr_tuples

        # Check for matching exponents and combine like terms
        new_arr = []
        for apar1, bpar1, epar1 in self.arr:
            found_like_term = False
            for i, (_, bpar2, epar2) in enumerate(new_arr):
                if numpy.allclose([bpar1, epar1], [bpar2, epar2]):
                    new_arr[i][0] += apar1
                    found_like_term = True

            if not found_like_term:
                new_arr.append([apar1, bpar1, epar1])
        self.arr = tuple(map(tuple, new_arr))

        self.arr_collid = collid

    def set_plog(self, plog_dct):
        """ Sets PLOG parameters

            :param plog_dct: dct describing PLOG parameters
            :type plog_dct: dct {pressure1: arr_tuples1, pressure2: ...}
        """

        self.check_plog(plog_dct)

        # If there are no PLOG parameters in the object, add them
        if self.plog is None:
            self.plog = {pressure: arr_tuples
                         for pressure, arr_tuples in plog_dct.items()
                         if pressure != 'high'}
        # If there are already PLOG parameters in the object, add dups
        else:
            if self.plog_dups is None:  # if no other dups
                self.plog_dups = [plog_dct]
            else:  # append new dups to other dups
                self.plog_dups.append(plog_dct)

    def set_cheb(self, cheb_dct):
        """ Sets Chebyshev parameters

            :param cheb_dct: dct describing Chebyshev parameters
            :type cheb_dct: dct {'tlim': (tmin, tmax), 'plim': (pmin, pmax),
                'alpha': matrix of Cheb coefficients}
        """

        self.check_cheb(cheb_dct)

        # Fill the one_atm_arr field if it's None (this may be redundant)
        if cheb_dct.get('one_atm_arr') is None:
            cheb_dct['one_atm_arr'] = None
        # If there are no Chebyshev parameters in the object, add them
        if self.cheb is None:
            self.cheb = cheb_dct
        # If there are already Chebyshev parameters in the object, add dups
        else:
            if self.cheb_dups is None:  # if no other dups
                self.cheb_dups = [cheb_dct]
            else:  # append new dups to other dups
                self.cheb_dups.append(cheb_dct)

    def set_troe(self, troe_dct):
        """ Sets Troe parameters

            :param troe_dct: dct describing Troe parameters
            :type troe_dct: dct {'highp_arr': highp_arr_tuples,
                'lowp_arr': lowp_arr_tuples,
                'troe_params': [alpha, T***, T*, T**],
                'collid': {spc1: eff1, spc2: ...}}
        """

        self.check_troe(troe_dct)

        # If there are no Troe parameters in the object, add them
        if self.troe is None:
            self.troe = troe_dct
        # If there are already Troe parameters in the object, add dups
        else:
            if self.troe_dups is None:  # if no other dups
                self.troe_dups = [troe_dct]
            else:  # append new dups to other dups
                self.troe_dups.append(troe_dct)

    def set_lind(self, lind_dct):
        """ Sets Lindemann parameters

            :param lind_dct: dct describing Lindemann parameters
            :type lind_dct: dct {'highp_arr': highp_arr_tuples,
                'lowp_arr': lowp_arr_tuples, 'collid': {spc1: eff1, spc2: ...}}
        """

        self.check_lind(lind_dct)

        # If there are no Lindemann parameters in the object, add them
        if self.lind is None:
            self.lind = lind_dct
        # If there are already Lindemann parameters in the object, add dups
        else:
            if self.lind_dups is None:  # if no other dups
                self.lind_dups = [lind_dct]
            else:  # append new dups to other dups
                self.lind_dups.append(lind_dct)

    def check_arr(self, arr_tuples):
        """ Ensures that Arrhenius parameters have correct form:
            ((A1, n1, Ea1), (A2, n2, Ea2), ...)

            :param arr_tuples: Arrhenius parameters
            :type arr_tuples: tuple(tuple(float))
        """

        assert isinstance(arr_tuples, (list, tuple)), (
            'arr_tuples should be a list-of-lists or tuple-of-tuples')
        for arr_tuple in arr_tuples:
            assert isinstance(arr_tuple, (list, tuple)), (
                'arr_tuple should be a list or tuple')
            assert len(arr_tuple) == 3, (
                'Length of each arr_tuple should be 3')
            assert all(isinstance(x, (int, float)) for x in arr_tuple), (
                'Each arr_tuple param should be a float or an int')

    def check_plog(self, plog_dct):
        """ Ensures that a PLOG dictionary has the correct form:
            {'high'/pressure: ((A1, n1, Ea1), (A2, n2, Ea2), ...),
             'high'/pressure: ...}

            :param plog_dct: dct describing PLOG parameters
            :type plog_dct: dct {pressure1: arr_tuples1, pressure2: ...}
        """

        assert isinstance(plog_dct, dict)
        for pressure, arr_tuples in plog_dct.items():
            assert pressure == 'high' or isinstance(pressure, (float, int)), (
                "PLOG pressures should be 'high', floats, or ints")
            self.check_arr(arr_tuples)

    def check_cheb(self, cheb_dct):
        """ Ensures that a Chebyshev dictionary is of the proper form:
            {'tlim': (tmin, tmax), 'plim': (pmin, pmax), 'alpha': alpha,
             'one_atm_arr': one_atm_arr}

            :param cheb_dct: dct describing Chebyshev parameters
            :type cheb_dct: dct {'tlim': (tmin, tmax), 'plim': (pmin, pmax),
                'alpha': matrix of Cheb coefficients}
        """

        tlim = cheb_dct.get('tlim')
        plim = cheb_dct.get('plim')
        alpha = cheb_dct.get('alpha')
        one_atm_arr = cheb_dct.get('one_atm_arr')

        if one_atm_arr is not None:  # the one atm rates are not required
            self.check_arr(one_atm_arr)
        assert isinstance(tlim, (list, tuple)), (
            f'Cheb tlim should be a list or tuple, not {type(tlim)}')
        assert isinstance(plim, (list, tuple)), (
            f'Cheb plim should be a list or tuple, not {type(plim)}')
        assert all(isinstance(x, (int, float)) for x in tlim), (
            'Cheb tlim entries should be floats or ints')
        assert all(isinstance(x, (int, float)) for x in plim), (
            'Cheb plim entries should be floats or ints')
        assert isinstance(alpha, (list, tuple, numpy.ndarray)), (
            'Cheb alpha should be lists, tuples, or Numpy arrays')
        for row in alpha:
            assert isinstance(row, (list, tuple, numpy.ndarray)), (
                'Cheb rows should be lists, tuples, or Numpy arrays')
            assert all(isinstance(x, (int, float)) for x in row), (
                'Cheb params should be floats or ints')

    def check_troe(self, troe_dct):
        """ Ensures that a Troe dictionary is of proper form:
            {'highp_arr': highp_arr, 'lowp_arr': lowp_arr,
             'troe_params': troe_params, 'collid': collid}

            :param troe_dct: dct describing Troe parameters
            :type troe_dct: dct {'highp_arr': highp_arr_tuples,
                'lowp_arr': lowp_arr_tuples,
                'troe_params': [alpha, T***, T*, T**],
                'collid': {spc1: eff1, spc2: ...}}
        """

        highp_arr = troe_dct['highp_arr']
        lowp_arr = troe_dct['lowp_arr']
        troe_params = troe_dct['troe_params']  # alpha, T***, T*, T**
        collid = troe_dct.get('collid')  # using get b/c it's optional

        self.check_arr(highp_arr)
        self.check_arr(lowp_arr)
        assert len(troe_params) in (3, 4), (
            f'Troe params should be 3 or 4 entries, not {len(troe_params)}')
        for troe_param in troe_params:
            if troe_param is not None:  # 4th entry can be None
                assert isinstance(troe_param, (float, int)), (
                    'Troe params should be floats or ints')
        if collid is not None:
            self.check_collid(collid)

    def check_lind(self, lind_dct):
        """ Ensures that a Lindemann dictionary is of proper form:
            {'highp_arr': highp_arr, 'lowp_arr': lowp_arr,
             'collid': collid}

            :param lind_dct: dct describing Lindemann parameters
            :type lind_dct: dct {'highp_arr': highp_arr_tuples,
                'lowp_arr': lowp_arr_tuples, 'collid': {spc1: eff1, spc2: ...}}
        """

        highp_arr = lind_dct['highp_arr']
        lowp_arr = lind_dct['lowp_arr']
        collid = lind_dct.get('collid')  # using get b/c it's optional

        self.check_arr(highp_arr)
        self.check_arr(lowp_arr)
        if collid is not None:
            self.check_collid(collid)

    def check_collid(self, collid):
        """ Ensures that a collider dictionary is of proper form:
            {spc1: eff1, spc2: ...}

            :param collid: collider efficiencies for various species
            :type collid: dct
        """

        assert isinstance(collid, dict), 'collid should be a dct'
        assert all(isinstance(spc, str) for spc in collid.keys()), (
            'collid keys (i.e., spcs) should be strings')
        assert all(isinstance(eff, (int, float)) for eff in collid.values()), (
            'collid values (i.e., efficiences) should be floats or ints')

    def get_existing_forms(self):
        """ Return a list of all the types of parameters contained
            in the object

            :return forms: all functional forms defined in the params
            :rtype: tuple
        """

        forms = ()
        if self.arr is not None:
            forms += ('arr',)
        if self.plog is not None:
            forms += ('plog',)
        if self.cheb is not None:
            forms += ('cheb',)
        if self.troe is not None:
            forms += ('troe',)
        if self.lind is not None:
            forms += ('lind',)

        return forms

    def combine_objects(self, other_params):
        """ Combines another RxnParams instance with the current one. This is
            used in the combining of duplicates when parsing mechanism files.

            :param other_params: another instance of RxnParams to be added
                to the current instance; should only have one functional form
            :type other_params: RxnParams object
        """

        # Get the forms to be added from the existing object
        form_to_be_added = other_params.get_existing_forms()
        assert len(form_to_be_added) == 1, (
            f'The object to be appended should'
            f' only have one functional form but has {form_to_be_added}.')

        # Combine the dictionaries depending on the form
        if 'arr' in form_to_be_added:
            new_arr = other_params.arr
            new_arr_dct = {'arr_tuples': new_arr}
            self.set_arr(new_arr_dct)
        elif 'plog' in form_to_be_added:
            new_plog = other_params.plog
            self.set_plog(new_plog)
        elif 'cheb' in form_to_be_added:
            new_cheb = other_params.cheb
            self.set_cheb(new_cheb)
        elif 'troe' in form_to_be_added:
            new_troe = other_params.troe
            self.set_troe(new_troe)
        elif 'lind' in form_to_be_added:
            new_lind = other_params.lind
            if self.lind_dups is None:
                self.lind_dups = [new_lind]
            else:
                self.lind_dups.append(new_lind)

    def check_for_dups(self):
        """ Checks for unsual cases of duplicates: either mixed forms (e.g.,
            PLOG and Arrhenius) or more than one of some form (e.g., two
            PLOGS).
            Note that this does not check for duplicate Arrhenius, as that info
            is readily contain by the arr attribute and is an acceptable case

            :return dups: whether or not any duplicates are present
            :rtype: Bool
            :return dup_counts: number of duplicates for each form (excluding
                Arrhenius). Does not include the single form; if a rxn has
                a single Chebyshev form, the Cheb dup_counts entry will be 0.
            :rtype: dct {'plog': num_plog_dups, 'cheb': ..., 'troe': ...,
                'lind': ...}
        """

        dups = False

        # The first sign of duplicates: if there is more than one form
        forms = self.get_existing_forms()
        if len(forms) > 1:
            dups = True

        # The second sign of duplicates: if there are any *_dups fields that
        # are not None
        any_dups = any((self.plog_dups, self.cheb_dups, self.troe_dups,
                        self.lind_dups))
        if any_dups:
            dups = True

        # Make the dup_counts dictionary
        dup_counts = {}
        if self.plog_dups is not None:
            dup_counts['plog'] = len(self.plog_dups)
        else:
            dup_counts['plog'] = 0
        if self.cheb_dups is not None:
            dup_counts['cheb'] = len(self.cheb_dups)
        else:
            dup_counts['cheb'] = 0
        if self.troe_dups is not None:
            dup_counts['troe'] = len(self.troe_dups)
        else:
            dup_counts['troe'] = 0
        if self.lind_dups is not None:
            dup_counts['lind'] = len(self.lind_dups)
        else:
            dup_counts['lind'] = 0

        return dups, dup_counts


def scale_factor(params, factor_A=None, factor_n=None, factor_EA=None):
    """ scale all rates by some factor
        A = A*factor_A
        n = n+factor_n
        EA = EA+factor_EA
    """
    def single_form(params, form, factor_A=None, factor_n=None, factor_EA=None):
        if form == 'arr':
            arr_tuples = params.arr
            arr_collid = params.arr_collid
            new_arr = fix_arr_tuples(
                arr_tuples, factor_A=factor_A, factor_n=factor_n, factor_EA=factor_EA)
            new_arr_dct = {'arr_tuples': new_arr, 'arr_collid': arr_collid}
            new_params = RxnParams(arr_dct=new_arr_dct)
        elif form == 'plog':
            plog_dct = params.plog
            new_plog_dct = {}
            for pressure, arr_tuples in plog_dct.items():
                new_arr = fix_arr_tuples(
                    arr_tuples, factor_A=factor_A, factor_n=factor_n, factor_EA=factor_EA)
                new_plog_dct[pressure] = new_arr
            new_params = RxnParams(plog_dct=new_plog_dct)
        elif form == 'lind':
            lind_dct = params.lind
            new_highp = fix_arr_tuples(
                lind_dct['highp_arr'], factor_A=factor_A, factor_n=factor_n, factor_EA=factor_EA)
            new_lowp = fix_arr_tuples(
                lind_dct['lowp_arr'], factor_A=factor_A, factor_n=factor_n, factor_EA=factor_EA)
            new_lind_dct = {'highp_arr': new_highp, 'lowp_arr': new_lowp,
                            'collid': lind_dct['collid']}
            new_params = RxnParams(lind_dct=new_lind_dct)
        elif form == 'troe':
            troe_dct = params.troe
            new_highp = fix_arr_tuples(
                troe_dct['highp_arr'], factor_A=factor_A, factor_n=factor_n, factor_EA=factor_EA)
            new_lowp = fix_arr_tuples(
                troe_dct['lowp_arr'], factor_A=factor_A, factor_n=factor_n, factor_EA=factor_EA)
            new_troe_dct = {'highp_arr': new_highp, 'lowp_arr': new_lowp,
                            'collid': troe_dct['collid'],
                            'troe_params': troe_dct['troe_params']}
            new_params = RxnParams(troe_dct=new_troe_dct)
        else:
            raise NotImplementedError(f'{form} form not working yet')

        return new_params

    def fix_arr_tuples(old_arr_tuples, factor_A=None, factor_n=None, factor_EA=None):
        new_arr_tuples = []
        for arr_tuple in old_arr_tuples:
            new_arr_tuple = list(arr_tuple)
            if factor_A:
                new_arr_tuple[0] *= factor_A
            if factor_n:
                new_arr_tuple[1] += factor_n
            if factor_EA:
                new_arr_tuple[2] += factor_EA
            new_arr_tuples.append(new_arr_tuple)

        return new_arr_tuples

    forms = params.get_existing_forms()
    for idx, form in enumerate(forms):
        if idx == 0:
            new_params = single_form(
                params, form, factor_A=factor_A, factor_n=factor_n, factor_EA=factor_EA)
        else:  # allows for more than one rate form; stupid, but general
            curr_params = single_form(
                params, form, factor_A=factor_A, factor_n=factor_n, factor_EA=factor_EA)
            new_params.combine_objects(curr_params)

    return new_params
