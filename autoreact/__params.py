import numpy
from autoreact.params import RxnParams

arr_dct = {'arr': [[1e12, 1.5, 50000], [1e12, 1.5, 50000]]}

params1 = RxnParams(arr_dct=arr_dct)

print(params1.arr)

plog_dct = {'high': [[1e12, 1.5, 50000], [1e12, 1.5, 50000]],
            1.0:    [[1e12, 1.5, 50000],],
            10.0:   [[1e12, 1.5, 50000], [1e12, 1.5, 50000]]
}


params2 = RxnParams(plog_dct=plog_dct)
print(params2.plog)

ref_alpha = numpy.array(
    [[1.86421309e+00, 4.20602838e-01, -5.74358452e-02, -5.45222311e-05],
     [7.61423648e+00, 7.51012552e-01, -9.23375204e-02, -8.25427040e-03],
     [-4.89211391e-01, 5.10360005e-01, -2.71105409e-02, -1.04075446e-02],
     [-3.93397030e-01, 2.67821927e-01, 1.58876205e-02, -6.32223880e-03],
     [-2.15290577e-01, 7.79192168e-02, 4.05605101e-02, 3.99721924e-03],
     [-8.40067735e-02, -2.24969601e-03, 2.50720909e-02, 5.36853083e-03]])

cheb_dct = {'alpha': ref_alpha,
            'tlim': (300, 2500),
            'plim': (0.01, 100),
            'one_atm_arr': [[1,0,0],]
}

params3 = RxnParams(cheb_dct=cheb_dct) 
