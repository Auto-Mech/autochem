"""Test autochem.rate."""

import pytest

from autochem import therm

NASA7 = """
! NASA polynomial coefficients for CH2O3
CH2O3(82)               H   2C   1O   3     G   100.000  5000.000  956.80      1
    1.18953474E+01 1.71602859E-03-1.47932622E-07 9.25919346E-11-1.38613705E-14    2
-7.81027323E+04-3.82545468E+01 2.98116506E+00 1.14388942E-02 2.77945768E-05    3
-4.94698411E-08 2.07998858E-11-7.51362668E+04 1.09457245E+01                   4
"""

MESSPF = """
Natural log of the partition function, its derivatives, entropy, and thermal capacity:
T, K    C5H7(487)z   C5H7(487)z   C5H7(487)z   C5H7(487)z   C5H7(487)z
               Z_0          Z_1          Z_2 S, cal/mol/K C, cal/mol/K
200        75.2225    0.0319355 -0.000107435        74.64      18.8322
300        78.0043    0.0246886 -4.72064e-05      82.9996       22.981
400        80.2842    0.0212802 -2.45403e-05      90.2987      28.0149
500        82.3078    0.0193399  -1.5554e-05      97.0643       32.692
600         84.172    0.0180099 -1.15225e-05      103.388       36.691
700        85.9194    0.0169743 -9.38185e-06      109.306      40.0754
800        87.5724    0.0161073 -8.04516e-06      114.851      42.9687
900        89.1446    0.0153523 -7.10332e-06      120.059      45.4678
1000       90.6456    0.0146793 -6.38525e-06      124.965      47.6397
1100       92.0826    0.0140705 -5.80925e-06      129.596      49.5327
1200       93.4614    0.0135142 -5.33095e-06      133.979      51.1852
1300       94.7869     0.013002 -4.92367e-06      138.134      52.6292
1400       96.0631    0.0125277 -4.57022e-06      142.082      53.8925
1500       97.2935    0.0120866 -4.25895e-06      145.838      54.9995
1600       98.4814    0.0116748 -3.98167e-06       149.42      55.9717
1700       99.6294    0.0112893 -3.73245e-06      152.839      56.8273
1800        100.74    0.0109275 -3.50685e-06      156.109      57.5825
1900       101.816    0.0105872 -3.30149e-06      159.241      58.2508
2000       102.858    0.0102666 -3.11369e-06      162.244      58.8441
2100        103.87   0.00996399 -2.94133e-06      165.128      59.3723
2200       104.851    0.0096779 -2.78262e-06      167.901      59.8441
2300       105.806   0.00940706  -2.6361e-06      170.571      60.2666
2400       106.733   0.00915031 -2.50052e-06      173.144      60.6461
2500       107.636   0.00890663  -2.3748e-06      175.627      60.9881
2600       108.515   0.00867505 -2.25802e-06      178.025      61.2971
2700       109.371   0.00845475 -2.14935e-06      180.343      61.5769
2800       110.206   0.00824494 -2.04807e-06      182.588      61.8311
2900       111.021   0.00804491 -1.95353e-06      184.761      62.0626
3000       111.816   0.00785403 -1.86517e-06      186.869      62.2739
298.2      77.9598    0.0247741 -4.78494e-05      82.8616      22.8932
"""


def generate_id(val):
    """Generate a unique ID for the test case."""
    if isinstance(val, str):
        return f"{val.strip()[:12]}"
    return val


@pytest.mark.parametrize(
    "name, spc_str0",
    [
        ("NASA7", NASA7),
    ],
    ids=generate_id,
)
def test__from_chemkin_string(name, spc_str0):
    # Dictionary serialization/deserialization
    spc = therm.from_chemkin_string(spc_str0)
    spc_ = therm.Species.model_validate(spc.model_dump())
    assert spc == spc_, f"\n   {spc}\n!= {spc_}"

    # Chemkin serialization/deserialization
    spc_str = therm.chemkin_string(spc)
    spc_ = therm.from_chemkin_string(spc_str)
    assert spc == spc_, f"\n   {spc}\n!= {spc_}"

    # Scale
    spc_times_2 = spc * 2
    spc_divided_by_2 = spc / 2

    # Plot
    therm.display(spc)
    therm.display(
        spc,
        label="original",
        others=[spc_times_2, spc_divided_by_2],
        others_labels=["doubled", "halved"],
    )


@pytest.mark.parametrize(
    "name, spc_str0",
    [
        ("MESSPF", MESSPF),
    ],
    ids=generate_id,
)
def test__fit(name, spc_str0):
    spc = therm.from_messpf_output_string(
        spc_str0,
        formula="C5H7",
        name="C5H7(487)z",
        Hf=86.239,
        Tf=0,
        units={"energy": "kcal"},
    )
    spc_fit = therm.fit(spc)
    therm.display(spc, label="data", others=[spc_fit], others_labels=["fit"])


if __name__ == "__main__":
    test__from_chemkin_string("NASA7", NASA7)
