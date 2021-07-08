"""  Functions to deal with vibrational frequencies
"""

import numpy
from phydat import phycon


def scale_frequencies_and_zpe(freqs, method, basis, scale_method='c3'):
    """ Scale frequencies according to some method
        obtain a corrected zpe
    """
    scaled_freqs = scale_frequencies(
        freqs,  method, basis, scale_method=scale_method)
    # Calculate the anharmonic zpe using the scaled freqs (anharm zpve)
    scaled_zpe = 0.0
    for freq, scfreq in zip(freqs, scaled_freqs):
        scaled_zpe += _anharm_zpve_from_scaling(freq, scfreq)
    scaled_zpe *= phycon.WAVEN2EH
    return scaled_freqs, scaled_zpe


def scale_frequencies(freqs, method, basis, scale_method='c3'):
    """ Scale frequencies according to some method
    """
    # Scale the frequencies
    if scale_method in SCALE_METHODS:
        scaled_freqs = SCALE_METHODS[scale_method](freqs, method, basis)
    else:
        scaled_freqs = freqs
    return scaled_freqs


def _anharm_zpve_from_scaling(freq, scaled_freq):
    """ Determine what the anharmonic ZPVE should be after scaling
    """
    return (freq / 2.0) - (1.0 / 8.0) * (scaled_freq - freq)


def rotor_scale_factor_from_harmonics(rt_freqs, rth_freqs, tors_freqs):
    """ scaling factor for rotor potentials to map them into harmonic
    """

    # Create a scaling factor for the frequencies
    # First sort tors frequencies in ascending order
    sort_tors_freqs = sorted(tors_freqs)

    # Keep only freqs whose RRHO freqs are above a threshold
    freq_thresh = 50.
    log_rt_freq = 0.0
    nfreq_remove = 0
    for freq in rt_freqs:
        if freq > freq_thresh:
            log_rt_freq += numpy.log(freq)
        else:
            nfreq_remove += 1

    log_freq = [numpy.log(freq) for freq in rth_freqs]
    log_freq = sum(log_freq)

    log_tors_freq = 0.0
    idx_remove = []
    for idx, freq in enumerate(sort_tors_freqs):
        if idx+1 > nfreq_remove:
            log_tors_freq += numpy.log(freq)
        else:
            idx_remove.append(tors_freqs.index(freq))

    # Generate the scaling factor
    factor = numpy.exp(log_rt_freq - log_freq - log_tors_freq)
    scale_factor = (idx_remove, factor)

    # generate the set of indices for torsions that are two be scales
    tau_factor = numpy.exp(log_rt_freq - log_freq)
    tau_factor_mode = tau_factor
    print(
        'TAU FACTOR {:4.6f} \t {:g} \t {:3.6f} {} '.format(
            tau_factor_mode, len(tors_freqs), factor,
            '-'.join([str(ridx) for ridx in idx_remove])))

    # Generate the set of indices for torsions that are two be scales
    scale_factor = (idx_remove, factor)

    return scale_factor


# Library of vibrational frequency scaling methods
M3_COEFFS_ANHARM = {
    ('b2plypd3', 'cc-pvtz'): (1.066, 0.008045, 0.33),
    ('wb97xd', '6-31g*'): (1.657244, 0.56000691, 0.029624)
    ('wb97xd', 'cc-pvtz'): (1.053471, 0.01186224, 0.26174883)
}

M3_COEFFS_HARM = {
    ('wb97xd', '6-31g*'): (0.91, -0.058, 0.001)
}


def _three_coeff_anharm_scaling(freqs, method, basis):
    """ Scales frequencies using factos with three coefficients
    """

    cf1, cf2, cf3 = M3_COEFFS_ANHARM.get((method, basis), (1.0, 0.0, 0.0))
    scaled_freqs = ()
    for freq in freqs:
        scale_factor = cf1 - (cf2 * freq**cf3)
        scaled_freqs += (freq * scale_factor,)

    return scaled_freqs


def _three_coeff_harm_scaling(freqs, method, basis):
    """ Scales frequencies using one factor, same factor applies to all frequencies
    """
    cf1, cf2, cf3 = M3_COEFFS_HARM.get((method, basis), (1.0, 0.0, 0.0))
    scaled_freqs = ()
    for freq in freqs:
        scale_factor = cf1 - (cf2 * freq**cf3)
        scaled_freqs += (freq * scale_factor,)
    return scaled_freqs


SCALE_METHODS = {
    'c3': _three_coeff_anharm_scaling,
    'c3_harm':  _three_coeff_harm_scaling
}
