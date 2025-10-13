"""testing autoreact functions."""

import numpy

import autoreact


def test_arr():
    """Tests various functions involving Arrhenius params."""
    # Test simple functions
    arr_dct = {"arr_tuples": [[1e12, 1.5, 50000], [1e12, 1.5, 50000]]}
    params1 = autoreact.params.RxnParams(arr_dct=arr_dct)
    assert params1.get_existing_forms() == ("arr",)
    assert numpy.allclose(params1.arr, [[2e12, 1.5, 50000]])
    dups, _ = params1.check_for_dups()
    assert dups is False  # should be False for Arrhenius


def test_plog():
    """Tests various functions involving PLOG params."""
    # Test simple functions
    plog_dct = {
        "high": [[1e12, 1.5, 50000], [1e12, 1.5, 50000]],
        1.0: [[1e12, 1.5, 50000]],
        10.0: [[1e12, 1.5, 50000], [1e12, 1.5, 50000]],
    }
    params2 = autoreact.params.RxnParams(plog_dct=plog_dct)
    assert params2.get_existing_forms() == ("plog",)
    for _, arr_tuples in params2.plog.items():
        for arr_tuple in arr_tuples:
            assert numpy.allclose([1e12, 1.5, 50000], arr_tuple)

    # Test duplicates
    params2.combine_objects(params2)  # combine params2 with itself
    assert len(params2.plog_dups) == 1
    dups, _ = params2.check_for_dups()
    assert dups


def test_cheb():
    """Tests various functions involving Chebyshev params."""
    # Test simple functions
    ref_alpha = numpy.array(
        [
            [1.86421309e00, 4.20602838e-01, -5.74358452e-02, -5.45222311e-05],
            [7.61423648e00, 7.51012552e-01, -9.23375204e-02, -8.25427040e-03],
            [-4.89211391e-01, 5.10360005e-01, -2.71105409e-02, -1.04075446e-02],
            [-3.93397030e-01, 2.67821927e-01, 1.58876205e-02, -6.32223880e-03],
            [-2.15290577e-01, 7.79192168e-02, 4.05605101e-02, 3.99721924e-03],
            [-8.40067735e-02, -2.24969601e-03, 2.50720909e-02, 5.36853083e-03],
        ]
    )
    cheb_dct = {
        "alpha": ref_alpha,
        "tlim": (300, 2500),
        "plim": (0.01, 100),
        "one_atm_arr": [[1, 0, 0]],
    }
    params3 = autoreact.params.RxnParams(cheb_dct=cheb_dct)
    assert params3.get_existing_forms() == ("cheb",)
    alpha = params3.cheb["alpha"]
    for idx, array in enumerate(alpha):
        assert numpy.allclose(array, ref_alpha[idx])

    # Test duplicates
    params3.combine_objects(params3)  # combine params3 with itself
    assert len(params3.cheb_dups) == 1
    dups, _ = params3.check_for_dups()
    assert dups


def test_troe():
    """Tests various functions involving Troe params."""
    troe_dct = {
        "highp_arr": [[1e12, 1.5, 50000]],
        "lowp_arr": [[1e12, 1.5, 50000]],
        "troe_params": [1.5, 8000, 100, 1000],
        "collid": {"AR": 1.4, "N2": 1.7},
    }
    params4 = autoreact.params.RxnParams(troe_dct=troe_dct)
    assert params4.get_existing_forms() == ("troe",)

    # Test duplicates
    params4.combine_objects(params4)  # combine params4 with itself
    assert len(params4.troe_dups) == 1
    dups, _ = params4.check_for_dups()
    assert dups


def test_lind():
    """Tests various functions involving Lindemann params."""
    lind_dct = {
        "highp_arr": [[1e12, 1.5, 50000]],
        "lowp_arr": [[1e12, 1.5, 50000]],
        "collid": {"AR": 1.4, "N2": 1.7},
    }
    params5 = autoreact.params.RxnParams(lind_dct=lind_dct)
    assert params5.get_existing_forms() == ("lind",)

    # Test duplicates
    params5.combine_objects(params5)  # combine params5 with itself
    assert len(params5.lind_dups) == 1
    dups, _ = params5.check_for_dups()
    assert dups


if __name__ == "__main__":
    test_arr()
    test_plog()
    test_cheb()
    test_troe()
    test_lind()
