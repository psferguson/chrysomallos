import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from dustmaps.sfd import SFDQuery

from chrysomallos.synthpop.artpop_source import MISTSersicSSPChrysomallos
from chrysomallos.utils import (
    mstar_from_absmag,
    rad_physical_to_sky,
    totmag_below_maglim,
)
from chrysomallos.utils.log import logger

__all__ = [
    "adopt_a_cat",
    "massage_the_cat",
]


def adopt_a_cat(
    wcs,
    bbox,
    age=10.0,
    feh=-2.0,
    stellar_mass=5.0e5,
    dist=2.0,
    r_scale=300.0,
    ellip=0,
    theta=0,
    n_sersic=1,
    m_v=None,
    mag_limit=36.0,
    mag_limit_band="LSST_g",
    random_seed=None,
    **kwargs,
):
    """Make a synthetic source catalog to inject into an image.

    Parameters
    ----------
    wcs : `wcs object`
        The wcs object associated with the image to inject into.
    bbox : `bbox object`
        The bounding box object associated with the image to inject into.
    age : `float`
        Age in Gyr
    feh : `float`
        [Fe/H]
    stellar_mass : `float`
        Total stellar mass in M_Sun
    dist : `float`
        Distance in Mpc
    r_scale : `float`
        Sersic scale radius in pc
    ellip: `float`
        ellipticity between 0:1. e=1-b/a
    theta: `float`
        position angle in deg
    n_sersic: `float`
        sersic index
    m_v: `float`
        Absolute V-band mag of dwarf to inject. Overrides "stellar_mass" if
        both are given as inputs. If not set, "stellar_mass" is a required input.
    mag_limit : `float`
        Faintest mag of stars to include
    mag_limit_band : `str`
        Band to apply mag_limit in
    random_seed: `int`
            if not None this sets the random seed to generate catalog with

    Returns
    -------
    cat : `Astropy Table`
        A catalog containing the simulated dwarf.
    """

    # use this random state for reproducibility
    np.random.seed(random_seed)

    # Note: we have assume rscale = r_half. I _think_ this is true for the
    #   projected half-light radius and Plummer scale radius...

    dist = dist * u.Mpc
    r_scale = r_scale * u.pc
    pixel_scale = wcs.getPixelScale().asArcseconds()
    # If the image is 2k x 2k, the dwarf will be centered at 1000, 1000
    xydim = 1999

    if m_v is not None:
        stellar_mass = mstar_from_absmag(m_v)

    # create the artpop stellar population
    ssp = MISTSersicSSPChrysomallos(
        log_age=np.log10(age * 1e9),
        feh=feh,
        total_mass=stellar_mass,
        distance=dist,
        r_eff=r_scale,
        ellip=ellip,
        theta=theta,
        n=n_sersic,
        phot_system="LSST",  # photometric system
        imf="kroupa",  # default imf
        xy_dim=xydim,  # half the size of an LSST patch
        pixel_scale=pixel_scale,  # pixel scale of input image
        mag_limit=mag_limit,
        mag_limit_band=mag_limit_band,
        add_remnants=False,  # don't include stellar remnants in the total mass calculation
    )

    x0 = bbox.beginX
    y0 = bbox.beginY
    xcoords = ssp.x + x0
    ycoords = ssp.y + y0
    radec_coords = wcs.pixelToSkyArray(xcoords, ycoords, degrees=True)

    # Add extinction
    # see https://github.com/rubin-dp0/delegate-contributions-dp02/blob/db7d0\
    # 6ba6203faa15732c5368ff6e52ea53c5796/MWhalo_density/Milky_Way_halo_\
    # density.ipynb#L64

    # set the A_lamba/E(B-V) values for the six ugrizy LSST filters
    band_a_ebv = np.array([4.812, 3.643, 2.699, 2.063, 1.578, 1.313])

    coords = SkyCoord(radec_coords[0], radec_coords[1], unit="deg", frame="icrs")

    sfd = SFDQuery()
    ebvvec = sfd(coords)

    u_ext = ssp.mags["LSST_u"] + ebvvec * band_a_ebv[0]
    g_ext = ssp.mags["LSST_g"] + ebvvec * band_a_ebv[1]
    r_ext = ssp.mags["LSST_r"] + ebvvec * band_a_ebv[2]
    i_ext = ssp.mags["LSST_i"] + ebvvec * band_a_ebv[3]
    z_ext = ssp.mags["LSST_z"] + ebvvec * band_a_ebv[4]

    # if band_for_injection == 'g':
    #     mag_for_injection = g_ext
    # elif band_for_injection == 'r':
    #     mag_for_injection = r_ext
    # elif band_for_injection == 'i':
    #     mag_for_injection = i_ext
    # else:
    #     mag_for_injection = i_ext

    dist_column = np.repeat(dist.value, len(ssp.mags))

    cat = Table(
        {
            "injection_id": np.arange(len(ssp.mags)),
            "ra": radec_coords[0],
            "dec": radec_coords[1],
            "source_type": ["DeltaFunction"] * len(ssp.mags),
            "dist": dist_column,
            "g_mag": g_ext,
            "r_mag": r_ext,
            "i_mag": i_ext,
            # 'mag': mag_for_injection,
        }
    )

    return cat


def massage_the_cat(
    cat_inp,
    mag_limit,
    band_for_injection,
    wcs,
    bbox,
    x_cen,
    y_cen,
    r_scale=300.0,
    dist=2.0,
    **kwargs,
):
    """Replace the total flux below some mag limit with the
         appropriate Sersic model.

    Parameters
    ----------
    cat : `Astropy Table`
        Catalog created by adopt_a_cat.py
    injection_maglim : `float`
        Magnitude limit below which the flux will be replaced by a Sersic model
    band_for_injection : `str`
        Band of the image that you're going to inject into
    x_cen : `float`
        X-coordinate to center the dwarf on (between 0-4000)
    y_cen : `float`
        Y-coordinate to center the dwarf on (between 0-4000)
    wcs : `wcs object`
        The wcs object associated with the image to inject into.
    bbox : `bbox object`
        The bounding box object associated with the image to inject into.
    r_scale : `float`
        Plummer scale radius in pc
    dist : `float`
        Distance in Mpc
    mag_limit : `float`
        Faintest mag of stars to include
    mag_limit_band : `str`
        Band to apply mag_limit in

    Returns
    -------
    cat : `Astropy Table`
        A catalog containing the simulated dwarf.
    """
    # check if catalog has any entries
    if len(cat_inp) < 1:
        raise RuntimeError("catalog must contain entries")
    # Calculate the total mag below the mag limit.
    # Reduce the original catalog to only stars above the magnitude limit.
    # Append the Sersic model for the remaining flux.

    band = band_for_injection + "_mag"

    mag_for_sersic = totmag_below_maglim(cat_inp[band], mag_limit)
    logger.debug(
        f"massage a cat mag limit {mag_limit}, sersic component {mag_for_sersic: 0.2f} mag"
    )
    # Replicate the magnitude column for the band you want to inject into:
    cat_inp["mag"] = cat_inp[band]

    cat = cat_inp[cat_inp["mag"] <= mag_limit]

    # Parameters for simulated dwarf:
    # r_scale = r_scale*u.pc
    # reff = (r_scale/1.68).to(u.kpc)
    # Get the angular size corresponding to the given radius, distance:
    radius = rad_physical_to_sky(r_scale, dist)
    pa = 0.0
    axis_ratio = 1.0
    # ra_sim = np.median(cat['ra'])
    # dec_sim = np.median(cat['dec'])

    # Shift the input coordinates to the desired xy coord:
    xycoords_inp = wcs.skyToPixelArray(cat["ra"], cat["dec"], degrees=True)
    # Using a fixed injection image size of 2000x2000 for now
    dwarfcen_x = 1000
    dwarfcen_y = 1000
    xxx = xycoords_inp[0] - dwarfcen_x + x_cen
    yyy = xycoords_inp[1] - dwarfcen_y + y_cen
    radec_coords_stars = wcs.pixelToSkyArray(xxx, yyy, degrees=True)
    cat["ra"] = radec_coords_stars[0]
    cat["dec"] = radec_coords_stars[1]

    # Append a single line for the "galaxy" that contains the unresolved flux:
    cat.add_row()

    x0 = bbox.beginX
    y0 = bbox.beginY
    xcoord = x_cen + x0
    ycoord = y_cen + y0
    radec_coords = wcs.pixelToSky(float(xcoord), float(ycoord))

    cat[-1]["ra"] = radec_coords[0].asDegrees()
    cat[-1]["dec"] = radec_coords[1].asDegrees()
    cat[-1]["mag"] = mag_for_sersic
    cat[-1]["source_type"] = "Sersic"

    # I think instead of including distance, I need to convert the radius
    # into a sky radius instead of physical. CHECK THIS!
    cat[-1]["dist"] = dist

    # We need Sersic index, position angle, axis ratio, and semimajor axis for
    # the galaxy model,
    #  so create columns for these (with all stars set to some default values)
    semimajor_all = 0.0 * cat["mag"]
    semimajor_all[-1] = radius
    sersic_n_all = 0.0 * cat["mag"] + 1.0
    sersic_n_all[-1] = 1.0
    pa_all = 0.0 * cat["mag"]
    axis_ratio_all = np.ones(len(cat["mag"]))

    cat.add_columns(
        [semimajor_all, sersic_n_all, pa_all, axis_ratio_all],
        names=["half_light_radius", "n_sersic", "pa", "axis_ratio"],
    )

    return cat
