import numpy as np

def totmag(mags):
    """Calculate the total magnitude of an input list of stars.

    Parameters
    ----------
    mags : `np.array`
        List of mags to sum.

    Returns
    -------
    mag_tot : `float`
        Total magnitude of all stars in input list.
    """

    # Take the first star in the list as a reference, and
    #   calculate fluxes relative to that reference star.
    mag_ref = mags[0]
    flux_vs_ref = 10.0**(-0.4*(mags-mag_ref))

    # Sum the fluxes of all the stars:
    flux_tot = np.sum(flux_vs_ref[1:])

    # Convert the summed flux to magnitude
    mag_tot = mag_ref - 2.5*np.log10(flux_tot)

    return mag_tot


def totmag_below_maglim(mags, maglim):
    """Calculate the total magnitude of an input list of stars.

    Parameters
    ----------
    mags : `np.array`
        List of mags to sum.
    maglim : `float`
        Sum only stars fainter than the maglim threshold.

    Returns
    -------
    mag_below_threshold : `float`
        Total magnitude of all stars with mag >= maglim from input list.
    """

    select_mags = (mags >= maglim)
    mag_below_threshold = totmag(mags[select_mags])

    return mag_below_threshold
