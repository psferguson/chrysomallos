import numpy as np

__all__ = [
    "totmag",
    "totmag_below_maglim",
    "fluxfrac_above_maglim",
    "mag_at_flux_percentile",
]


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
    flux_tot = np.sum(flux_vs_ref)

    # Convert the summed flux to magnitude
    mag_tot = mag_ref - 2.5*np.log10(flux_tot)

    return mag_tot


def totmag_below_maglim(mags, maglim):
    """Calculate the total magnitude of an input list of stars
         below some limiting mag.

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


def fluxfrac_above_maglim(mags, maglim):
    """Calculate the fraction of the flux from an input list of stars
         above some limiting mag.

    Parameters
    ----------
    mags : `np.array`
        List of magnitudes of the stars.
    maglim : `float`
        Sum the flux of stars brighter than the maglim threshold.

    Returns
    -------
    frac_above_maglim : `float`
        Fraction of the total flux of input star list that comes from
        mag >= maglim.
    """

    sorted_mags = np.sort(mags)
    mag_ref = sorted_mags[0]
    flux_vs_ref = 10.0**(-0.4*(sorted_mags-mag_ref))
    flux_tot = np.sum(flux_vs_ref)
    select_mags = (sorted_mags >= maglim)
    frac_above_maglim = 1.0 - np.sum(flux_vs_ref[select_mags])/flux_tot

    return frac_above_maglim


def mag_at_flux_percentile(mags, pct):
    """Calculate the magnitude at which the fraction of the
         cumulative flux is equal to the given percentile.
         (e.g., above what mag is 90% of the flux contained?)

    Parameters
    ----------
    mags : `np.array`
        List of magnitudes of the stars.
    pct : `float`
        Percentage of flux for which you want the corresponding mag.

    Returns
    -------
    mag_at_pct : `float`
        Magnitude above which the requested fraction of the total flux
        is contained.
    """

    sorted_mags = np.sort(mags)
    mag_ref = sorted_mags[0]
    flux_vs_ref = 10.0**(-0.4*(sorted_mags-mag_ref))
    flux_tot = np.sum(flux_vs_ref)
    cumflux_frac = np.cumsum(flux_vs_ref)/flux_tot
    mag_at_pct = sorted_mags[np.argmin(np.abs(cumflux_frac-pct))]

    return mag_at_pct
