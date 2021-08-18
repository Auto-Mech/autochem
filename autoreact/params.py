class RxnParams:
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
        assert len(arr_params) == 3, ('Length of arr_params should be 3 but'
            ' is {len(arr_params)}')

        # Add the arr_params to the RxnParams object
        if self.arr is None:
            self.arr = [arr_params]
        else:
            self.arr.append(arr_params)

        # If a second set of Arrhenius parameters was input
        if arr_params2 is not None:
            assert isinstance(arr_params2, list), 'arr_params2 should be a list'
            assert len(arr_params2) == 3, ('Length of arr_params2 should be 3'
                ' but is {len(arr_params2)}')
            self.arr.append(arr_params2)


    def set_plog():
        pass

    def set_cheb():
        pass

    def set_troe():
        pass            

