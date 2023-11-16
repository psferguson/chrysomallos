import numpy as np
import astropy.units as u

__all__ = [
    "totmag",
    "totmag_below_maglim",
    "fluxfrac_above_maglim",
    "mag_at_flux_percentile",
    "rad_physical_to_sky",
    "get_flux_in_annulus",
    "sb_rh_to_mv",
    "sb_mv_to_rh",
    "mstar_from_absmag"
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
    if select_mags.sum() > 1:
        mag_below_threshold = totmag(mags[select_mags])
    else:
        return 40
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


def rad_physical_to_sky(radius, distance):
    """Convert a physical size at a given distance to angular size on the sky.

    Parameters
    ----------
    radius : `float`
        Radius (size) you wish to convert, in pc
    distance : `float`
        Distance in Mpc

    Returns
    -------
    rad_arcsec : `float`
        Radius converted to arseconds.
    """
    radius = radius*u.pc
    distance = distance*u.Mpc
    angle = radius/(distance.to(u.pc))*u.rad
    return angle.to(u.arcsec).value


def get_flux_in_annulus(image, xpos, ypos, r_inner, r_outer):
    """
    Extract the flux within an annulus centered on a given position from an image.

    Parameters
    ----------
    image: ExposureF
        Image to measure.
    xpos: float
        X pixel position to center on.
    ypos: float
        Y pixel position to center on.
    r_inner: float
        Inner radius of annulus, in arcseconds
    r_outer: float
        Outer radius of annulus, in arcseconds

    Returns
    -------
    totflux
    """
    pix_scale = image.getWcs().getPixelScale().asArcseconds()
    r_inner_pix = r_inner / pix_scale
    r_outer_pix = r_outer / pix_scale
    xv = np.arange(0, image.getWidth(), 1)
    yv = np.arange(0, image.getHeight(), 1)
    xx, yy = np.meshgrid(xv, yv)
    rad = np.sqrt((xx-xpos)**2 + (yy-ypos)**2)
    picksel = (rad > r_inner_pix) & (rad < r_outer_pix)
    totflux = np.sum(image.image.array[picksel])
    return totflux


def absmag_from_mstar(mstars, m_to_l=1.6):
    """
    From an input satellite stellar mass, calculate the absolute magnitude.
    Assuming M*/LV = 1.6 for dSphs (Woo+2008).

    L_V/L_Sun = 10^[(M_V,Sun-M_V)/2.5]

    Parameters
    ----------
    m_v : `float`
        Luminosity (V-band absolute magnitude) of the satellite
    m_to_l : `float` (default: 1.6)
        Stellar mass to (V-band) light for a dSph.

    Returns
    -------
    mstars : `float`
        Stellar mass of the satellite in solar masses

    """

    mv_sun = 4.83

    lv = mstars / m_to_l
    m_v = mv_sun - 2.5*np.log10(lv)

    return m_v


def mstar_from_absmag(m_v, m_to_l=1.6):
    """
    From an input satellite M_V, calculate the luminosity in solar units.
    Assuming M*/LV = 1.6 for dSphs (Woo+2008), infer the stellar mass.

    L_V/L_Sun = 10^[(M_V,Sun-M_V)/2.5]

    Parameters
    ----------
    m_v : `float`
        Luminosity (V-band absolute magnitude) of the satellite
    m_to_l : `float` (default: 1.6)
        Stellar mass to (V-band) light for a dSph.

    Returns
    -------
    mstars : `float`
        Stellar mass of the satellite in solar masses

    """

    mv_sun = 4.83

    lv = 10.0**((mv_sun-m_v)/2.5)
    mstars = m_to_l * lv

    return mstars


def sb_rh_to_mv(sb, rh, distance):
    """
    From an input surface brightness and half-light radius, calculate
    the absolute magnitude. Assumes a circular dwarf (i.e., radius, not a).

    Parameters
    ----------
    sb : `float`, mag/arcsec**2
        Surface brightness within r_half (in mag/arcsec**2) of the satellite
    rh : `float`, pc
        Half-light radius (in pc) of the satellite
    distance : `float`, pc
        Distance to the satellite (in pc)

    Returns
    -------
    M_v : `float`
        Absolute luminosity (V-band absolute magnitude) of the satellite

    """

    r_over_d_radians = rh/distance
    r_over_d_arcsec = np.rad2deg(r_over_d_radians)*3600.0
    area_arcsec = np.pi * (r_over_d_arcsec**2)
    mv = sb - 2.5*np.log10(area_arcsec)
    M_V = mv - 5.0*np.log10(distance) + 5.0
    return M_V


def sb_mv_to_rh(sb, M_v, distance):
    """
    From an input surface brightness and half-light radius, calculate
    the absolute magnitude. Assumes a circular dwarf (i.e., radius, not a).

    Parameters
    ----------
    sb : `float`, mag/arcsec**2
        Surface brightness within r_half (in mag/arcsec**2) of the satellite
    M_v : `float`
        Absolute luminosity (V-band absolute magnitude) of the satellite
    distance : `float`, pc
        Distance to the satellite (in pc)

    Returns
    -------
    rh : `float`, pc
        Half-light radius (in pc) of the satellite

    """
    mv = M_v + 5.0*np.log10(distance) - 5.0
    rh_over_distance_arcsec = np.sqrt((10.0**((sb-mv)/2.5))/np.pi)
    rh_over_distance_deg = rh_over_distance_arcsec/3600.0
    rh_over_distance_rad = np.deg2rad(rh_over_distance_deg)
    rh = rh_over_distance_rad * distance
    return rh
