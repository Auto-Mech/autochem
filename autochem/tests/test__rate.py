"""Test autochem.rate."""

import numpy
import pytest

from autochem import rate

# Elementary
SIMPLE = {
    "units": {"energy": "kcal"},
    "chemkin": """
H(4)+HO2(10)=H2(2)+O2(3)                            2.945e+06 2.087     -1.455   
""",
}

# Three-body (now same type as Elementary in Cantera)
THREEBODY = {
    "units": {"energy": "cal"},
    "chemkin": """
H2(2)+M=H(4)+H(4)+M                                 4.580e+19 -1.400    104.40
H2(2)/2.55/ N2/1.01/ CO(1)/1.95/ CO2(12)/3.83/ H2O(7)/12.02/ 
""",
}

# Falloff
FALLOFF = {
    "units": {"energy": "kcal"},
    "chemkin": """
CO(1)+O(5)(+M)=CO2(12)(+M)                          1.880e+11 0.000     2.430    
CO2(12)/3.80/ CO(1)/1.90/ H2O(7)/12.00/ 
    LOW/ 1.400e+21 -2.100    5.500    /
""",
}

# Falloff (Troe)
FALLOFF_TROE = {
    "units": {"energy": "cal"},
    "chemkin": """
C2H4(+M)=H2+H2CC(+M)     8.000E+12     0.440    88770
    LOW  /               7.000E+50    -9.310    99860  /
    TROE /   7.345E-01   1.800E+02   1.035E+03   5.417E+03 /
    H2/2.000/   H2O/6.000/   CH4/2.000/   CO/1.500/  CO2/2.000/  C2H6/3.000/  AR/0.700/
""",
}

# Activated
ACTIVATED = {
    "units": {"energy": "cal"},
    "chemkin": """
CH3+CH3(+M)=H + C2H5(+M) 4.989E12 0.099 10600.0 ! Stewart
HIGH/ 3.80E-7 4.838 7710. / ! Chemically activated reaction
""",
}

# Activated (SRI)
ACTIVATED_SRI = {
    "units": {"energy": "cal"},
    "chemkin": """
CH3+CH3(+M)=H + C2H5(+M) 4.989E12 0.099 10600.0 ! Stewart
HIGH/ 3.80E-7 4.838 7710. / ! Chemically activated reaction
SRI / 1.641 4334 2725 / ! SRI pressure dependence
""",
}

# Plog
PLOG = {
    "units": {"energy": "cal"},
    "chemkin": """
C2H4+OH=PC2H4OH          2.560E+36    -7.752     6946   ! Arrhenius parameters at 1 atm
    PLOG /   1.000E-02   1.740E+43   -10.460     7699 /
    PLOG /   2.500E-02   3.250E+37    -8.629     5215 /
    PLOG /   1.000E-01   1.840E+35    -7.750     4909 /
    PLOG /   1.000E+00   2.560E+36    -7.752     6946 /
    PLOG /   1.000E+01   3.700E+33    -6.573     7606 /
    PLOG /   1.000E+02   1.120E+26    -4.101     5757 /
""",
}

# Cheb
CHEB = {
    "units": {"energy": "kcal"},
    "chemkin": """
CHO2(38)(+M)=H(4)+CO2(12)(+M)                       1.000e+00 0.000     0.000    
    TCHEB/ 300.000   2200.000 /
    PCHEB/ 0.010     98.702   /
    CHEB/ 6 4/
    CHEB/ 4.859e+00    9.247e-01    -3.655e-02   2.971e-15   /
    CHEB/ -3.682e-02   4.836e-02    2.294e-02    4.705e-17   /
    CHEB/ -3.061e-15   -4.107e-16   1.623e-17    -1.319e-30  /
    CHEB/ 3.682e-02    -4.836e-02   -2.294e-02   -4.705e-17  /
    CHEB/ -4.859e+00   -9.247e-01   3.655e-02    -2.971e-15  /
    CHEB/ 3.682e-02    -4.836e-02   -2.294e-02   -4.705e-17  /
""",
}

MESS1 = {
    "data": r"""
W1->W2

P\T           500       600       700       800       900     1e+03   1.1e+03   1.2e+03   1.3e+03   1.4e+03   1.5e+03   1.6e+03   1.7e+03   1.8e+03   1.9e+03     2e+03
0.1       0.00978      3.87       257  5.19e+03       ***  1.98e+05  5.83e+05       ***       ***       ***       ***       ***       ***       ***       ***       ***
1         0.00988      4.07       301  7.29e+03  7.46e+04  4.72e+05  1.77e+06  4.72e+06       ***       ***       ***       ***       ***       ***       ***       ***
10        0.00989       4.1       310  8.01e+03  9.81e+04  7.07e+05   3.3e+06   1.1e+07  2.78e+07       ***       ***       ***       ***       ***       ***       ***
100       0.00989       4.1       312  8.13e+03  1.03e+05  7.87e+05  4.09e+06  1.58e+07  4.74e+07  1.17e+08  2.43e+08  4.45e+08  7.33e+08  1.11e+09       ***       ***
1e+03     0.00989       4.1       312  8.14e+03  1.04e+05     8e+05  4.26e+06  1.71e+07  5.52e+07  1.49e+08  3.47e+08  7.15e+08  1.33e+09  2.25e+09  3.54e+09  5.23e+09
1e+04     0.00989       4.1       312  8.14e+03  1.04e+05  8.01e+05  4.28e+06  1.73e+07  5.64e+07  1.55e+08  3.73e+08  8.02e+08  1.57e+09  2.83e+09  4.78e+09   7.6e+09
1e+05     0.00989       4.1       312  8.14e+03  1.04e+05  8.01e+05  4.28e+06  1.73e+07  5.65e+07  1.56e+08  3.77e+08  8.14e+08  1.61e+09  2.94e+09  5.04e+09  8.19e+09
1e+06     0.00989       4.1       312  8.14e+03  1.04e+05  8.02e+05  4.28e+06  1.73e+07  5.66e+07  1.56e+08  3.77e+08  8.15e+08  1.61e+09  2.95e+09  5.07e+09  8.26e+09
1e+07     0.00989       4.1       312  8.14e+03  1.04e+05  8.02e+05  4.28e+06  1.73e+07  5.66e+07  1.56e+08  3.77e+08  8.16e+08  1.61e+09  2.95e+09  5.08e+09  8.27e+09
O-O       0.00991      4.11       312  8.15e+03  1.04e+05  8.02e+05  4.28e+06  1.73e+07  5.66e+07  1.56e+08  3.77e+08  8.16e+08  1.61e+09  2.95e+09  5.08e+09  8.28e+09
""",
    "fit": """
C5H9(553) = C5H9(536)              1.000      0.000      0.000   # pes.subpes.channel  1.1.1
    PLOG  /    0.01000  7.030E+74     -19.80      52500/
    PLOG  /     0.1000  5.530E+73     -19.00      55000/
    PLOG  /      1.000  2.010E+63     -15.50      52500/
    PLOG  /      5.000  6.020E+50     -11.50      48500/
    PLOG  /      10.00  1.510E+46     -10.00      47300/
    PLOG  /      20.00  6.020E+44     -9.500      47300/
""",
}

MESS2 = {
    "data": r"""
W2->W14

P\T           500       600       700       800       900     1e+03   1.1e+03   1.2e+03   1.3e+03   1.4e+03   1.5e+03   1.6e+03   1.7e+03   1.8e+03   1.9e+03     2e+03
0.1           ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***
1             ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***
10            ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***
100           ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***
1e+03         ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***
1e+04    1.69e+04  3.77e+05  3.49e+06       ***  6.86e+07  1.96e+08  4.63e+08  9.47e+08  1.74e+09  2.92e+09  4.59e+09  6.78e+09  9.56e+09   1.3e+10   1.7e+10  2.16e+10
1e+05    1.74e+04   3.8e+05  3.48e+06  1.85e+07  6.83e+07  1.95e+08  4.62e+08   9.5e+08  1.75e+09  2.97e+09  4.71e+09  7.04e+09  1.01e+10  1.38e+10  1.84e+10  2.38e+10
1e+06    1.75e+04  3.81e+05  3.48e+06  1.85e+07  6.82e+07  1.95e+08  4.62e+08   9.5e+08  1.76e+09  2.98e+09  4.72e+09  7.07e+09  1.01e+10   1.4e+10  1.86e+10  2.41e+10
1e+07    1.74e+04  3.81e+05  3.49e+06  1.85e+07  6.82e+07  1.95e+08  4.61e+08   9.5e+08  1.76e+09  2.98e+09  4.72e+09  7.08e+09  1.01e+10   1.4e+10  1.86e+10  2.42e+10
O-O      1.75e+04  3.81e+05  3.49e+06  1.85e+07  6.83e+07  1.95e+08  4.62e+08  9.51e+08  1.76e+09  2.98e+09  4.72e+09  7.08e+09  1.01e+10   1.4e+10  1.86e+10  2.42e+10
""",
}


@pytest.mark.parametrize(
    "name, data, check_roundtrip",
    [
        ("SIMPLE", SIMPLE, True),
        ("THREEBODY", THREEBODY, True),
        ("FALLOFF", FALLOFF, True),
        ("FALLOFF_TROE", FALLOFF_TROE, True),
        ("ACTIVATED", ACTIVATED, True),
        ("ACTIVATED_SRI", ACTIVATED_SRI, True),
        ("PLOG", PLOG, True),
        ("CHEB", CHEB, False),
    ],
)
def test__from_chemkin_string(name, data, check_roundtrip: bool):
    units = data.get("units")
    rxn_str0 = data.get("chemkin")

    # Dictionary serialization/deserialization
    rxn = rate.from_chemkin_string(rxn_str0, units=units)
    rxn_ = rate.Reaction.model_validate(rxn.model_dump())
    if check_roundtrip:
        assert rxn == rxn_, f"\n   {rxn}\n!= {rxn_}"

    # Chemkin serialization/deserialization
    rxn_str = rate.chemkin_string(rxn)
    rxn_ = rate.from_chemkin_string(rxn_str)
    if check_roundtrip:
        assert rxn == rxn_, f"\n   {rxn}\n!= {rxn_}"

    # Scale
    rxn_times_2 = rxn * 2
    rxn_divided_by_2 = rxn / 2

    # Plot
    rate.display(rxn)
    rate.display(
        [rxn, rxn_times_2, rxn_divided_by_2],
        label=["original", "doubled", "halved"],
    )

    # Evaluate
    T0 = 500
    T1 = [500, 600, 700, 800]
    P0 = 1.0
    P1 = [0.1, 1.0, 10.0]
    kT0P0 = rxn.rate(T0, P0)
    assert numpy.shape(kT0P0) == (), kT0P0
    kT1P0 = rxn.rate(T1, P0)
    assert numpy.shape(kT1P0) == (4,), kT1P0
    kT0P1 = rxn.rate(T0, P1)
    assert numpy.shape(kT0P1) == (3,), kT0P1
    kT1P1 = rxn.rate(T1, P1)
    assert numpy.shape(kT1P1) == (4, 3), kT1P1


@pytest.mark.parametrize(
    "name, data",
    [
        ("MESS1", MESS1),
        ("MESS2", MESS2),
    ],
)
def test__from_mess_channel_output(name, data):
    mess_chan_out = data.get("data")
    rxn_fit_str = data.get("fit")
    rxn = rate.from_mess_channel_output(mess_chan_out)

    # Scale
    rxn_times_2 = rxn * 2
    rxn_divided_by_2 = rxn / 2

    # Plot
    rate.display(rxn)
    rate.display(rxn, label="k")
    rate.display(
        [rxn, rxn_times_2, rxn_divided_by_2],
        label=["original", "doubled", "halved"],
    )

    if rxn_fit_str is not None:
        rxn_fit = rate.from_chemkin_string(rxn_fit_str)
        rate.display(
            [rxn, rxn_times_2, rxn_divided_by_2, rxn_fit],
            label=["original", "doubled", "halved", "fit"],
        )

    # Evaluate
    T0 = 500
    T1 = [500, 600, 700, 800]
    P0 = 1.0
    P1 = [0.1, 1.0, 10.0]
    kT0P0 = rxn.rate(T0, P0)
    assert numpy.shape(kT0P0) == (), kT0P0
    kT1P0 = rxn.rate(T1, P0)
    assert numpy.shape(kT1P0) == (4,), kT1P0
    kT0P1 = rxn.rate(T0, P1)
    assert numpy.shape(kT0P1) == (3,), kT0P1
    kT1P1 = rxn.rate(T1, P1)
    assert numpy.shape(kT1P1) == (4, 3), kT1P1


@pytest.mark.parametrize(
    "name, data",
    [
        ("MESS1", MESS1),
        ("MESS2", MESS2),
    ],
)
def test__fit_high(name, data):
    mess_chan_out = data.get("data")
    rxn = rate.from_mess_channel_output(mess_chan_out)
    rxn_fit = rate.fit_high(rxn)
    rate.display([rxn, rxn_fit], label=["rate", "fit"])


def test__fit_plog():
    # S(1210)r0 = OH(4) + S(1288)rs0 (1,2-epoxy)
    e12_mess_out_fake = r"""
    W2->W15

    P\T           500       600       700       800       900     1e+03   1.1e+03   1.2e+03   1.3e+03   1.4e+03   1.5e+03   1.6e+03   1.7e+03   1.8e+03   1.9e+03     2e+03
    0.1         0.692      16.9      98.5       300  1.11e+03  1.64e+04  1.72e+05  1.03e+06  4.42e+06  1.26e+07  2.93e+07  6.25e+07  1.14e+08  1.82e+08  2.81e+08  4.08e+08
    1            1.91      56.1       517  3.13e+03  2.25e+04  2.17e+05  1.38e+06  5.93e+06     2e+07  4.84e+07     1e+08  1.92e+08  3.24e+08  4.88e+08  7.14e+08  9.82e+08
    10           49.3  1.44e+03  1.81e+04  1.45e+05   7.6e+05  3.26e+06   1.1e+07  3.07e+07  7.62e+07  1.54e+08  2.81e+08  4.85e+08  7.71e+08  1.11e+09  1.58e+09   2.1e+09
    100      5.42e+03  1.11e+05  9.82e+05  5.27e+06  1.93e+07  5.42e+07  1.22e+08  2.31e+08  3.99e+08  6.17e+08  9.13e+08  1.31e+09  1.84e+09  2.45e+09  3.24e+09  4.13e+09
    O-O      1.75e+04  3.81e+05  3.49e+06  1.85e+07  6.83e+07  1.95e+08  4.62e+08  9.51e+08  1.76e+09  2.98e+09  4.72e+09  7.08e+09  1.01e+10   1.4e+10  1.86e+10  2.42e+10
    """
    e12_mess_out_real = r"""
    W2->P3

    P\T           500       600       700       800       900     1e+03   1.1e+03   1.2e+03   1.3e+03   1.4e+03   1.5e+03   1.6e+03   1.7e+03   1.8e+03   1.9e+03     2e+03
    0.1      1.46e+04  2.13e+05  1.02e+06  2.55e+06  4.77e+06  8.66e+06  1.61e+07  2.92e+07  4.99e+07  7.26e+07  1.04e+08   1.3e+08  1.64e+08  1.91e+08  2.18e+08  2.16e+08
    1         1.7e+04   3.3e+05   2.3e+06  7.85e+06   1.7e+07  2.94e+07  4.64e+07  7.02e+07  1.04e+08  1.39e+08  1.87e+08  2.21e+08  2.65e+08  2.99e+08   3.3e+08  3.16e+08
    10       1.74e+04  3.71e+05  3.19e+06  1.46e+07  4.13e+07   8.3e+07  1.33e+08  1.85e+08  2.37e+08  2.81e+08  3.32e+08  3.59e+08  3.99e+08  4.31e+08  4.58e+08  4.26e+08
    100      1.21e+04  2.69e+05  2.47e+06  1.26e+07  4.29e+07  1.05e+08     2e+08   3.1e+08  4.11e+08  4.95e+08  5.66e+08  5.93e+08  6.22e+08  6.44e+08  6.54e+08  6.08e+08
    O-O           ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***       ***
    """
    e12_rxn_fake = rate.from_mess_channel_output(e12_mess_out_fake)
    e12_rxn_real = rate.from_mess_channel_output(e12_mess_out_real)
    e12_rxn_data = e12_rxn_fake + e12_rxn_real
    e12_rxn_fit = rate.fit_plog(e12_rxn_data)
    rate.display([e12_rxn_fake, e12_rxn_real, e12_rxn_data, e12_rxn_fit])


@pytest.mark.parametrize(
    "name, data, exp_dct, count, factor",
    [
        (
            "PLOG",
            PLOG,
            {
                "C2H4": ["C2H4a", "C2H4b", "C2H4c"],
                "OH": ["OHa", "OHb", "OHc", "OHd"],
                "PC2H4OH": ["PC2H4OHa", "PC2H4OHb"],
            },
            24,
            1 / 2.0,
        ),
    ],
)
def test__expand_lumped(name, data, exp_dct, count: int, factor: float):
    units = data.get("units")
    rxn_str0 = data.get("chemkin")

    # Read
    rxn = rate.from_chemkin_string(rxn_str0, units=units)
    rxns = rate.expand_lumped(rxn, exp_dct)
    assert len(rxns) == count, f"{len(rxns)} != {count}"

    assert all(
        rxn.rate * factor == r.rate for r in rxns
    ), f"{rxn.rate} * {factor} != {rxns[0].rate}"


if __name__ == "__main__":
    # test__from_chemkin_string("FALLOFF", FALLOFF)
    # test__from_chemkin_string("THREEBODY", THREEBODY)
    # test__from_chemkin_string("PLOG", PLOG)
    # test__from_chemkin_string("CHEB", CHEB, False)
    test__expand_lumped(
        "PLOG",
        PLOG,
        {
            "C2H4": ["C2H4a", "C2H4b", "C2H4c"],
            "OH": ["OHa", "OHb", "OHc", "OHd"],
            "PC2H4OH": ["PC2H4OHa", "PC2H4OHb"],
        },
        24,
        1 / 2.0,
    )
