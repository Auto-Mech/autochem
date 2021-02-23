""" tunneling
"""

ENES = (1.0, 2.0)
ALPHA = 1.0  # value that goes with the values calculated below
RXN_FREQ = 1.0

# alpha
CFC_DCT = {}
QFC_DCT = {}
FREQS = []
RIDX = 0


def test__():
    """ test automol.reac.tunnel.transmission_coefficients
        test automol.reac.tunnel.actions
    """

    ref_alpha1 = 1.0
    ref_alpha2 = ALPHA

    alpha1 = automol.reac.tunnel.alpha(freqs, cfc_dct)
    alpha2 = automol.reac.tunnel.alpha(freqs, cfc_dct, qfc_dct=QFC_DCT)
    print(alpha1)
    print(alpha2)
    # assert numpy.isclose(alpha1, ref_alpha1)
    # assert numpy.isclose(alpha2, ref_alpha2)

    ref_kes = ()
    ref_ses = ()

    kes = automol.reac.tunnel.transmission_coefficient(enes, alpha, rxn_freq)
    ses = automol.reac.tunnel.action(enes, alpha, rxn_freq)
    # assert numpy.allclose(kes, ref_kes)
    # assert numpy.allclose(ses, ref_ses)
    print(kes)
    print(ses)


if __name__ == '__main__':
    test__()
