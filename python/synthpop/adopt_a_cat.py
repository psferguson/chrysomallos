import artpop
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
import numpy as np
import pandas as pd
from utils import totmag, totmag_below_maglim, fluxfrac_above_maglim, mag_at_flux_percentile

def adopt_a_cat(wcs, bbox, xcen, ycen, band_for_injection,
                age=10.0, feh=-2.0, mass=5.0e5, dist=2.0, r_scale=300.0,
                mag_limit=36.0, mag_limit_band='LSST_g',
                ):
    """Make a synthetic source catalog to inject into an image.

    Parameters
    ----------
    wcs : `wcs object`
        The wcs object associated with the image to inject into.
    bbox : `bbox object`
        The bounding box object associated with the image to inject into.
    xcen : `float`
        X-coordinate to center the dwarf on (between 0-4000)
    ycen : `float`
        Y-coordinate to center the dwarf on (between 0-4000)
    band_for_injection : `str`
        Band of image that the stars will be injected into.
    age : `float`
        Age in Gyr
    feh : `float`
        [Fe/H]
    mass : `float`
        Total stellar mass in M_Sun
    dist : `float`
        Distance in Mpc
    r_scale : `float`
        Plummer scale radius in pc
    mag_limit : `float`
        Faintest mag of stars to include
    mag_limit_band : `str`
        Band to apply mag_limit in

    Returns
    -------
    cat : `Astropy Table`
        A catalog containing the simulated dwarf.
    """

    # use this random state for reproducibility
    rand = np.random.RandomState()

    # Note: we have assume rscale = r_half. I _think_ this is true for the
    #   projected half-light radius and Plummer scale radius...

    dist = dist*u.Mpc
    r_scale = r_scale*u.pc

    # create the artpop stellar population
    ssp = artpop.MISTPlummerSSP(
        log_age=np.log10(age*1e9),
        feh=feh,
        total_mass=mass,
        distance=dist,
        scale_radius=r_scale,
        phot_system='LSST', # photometric system
        imf='kroupa', # default imf
        random_state=rand, # random state (can be set for reproducibility)
        xy_dim=1999, # half the size of an LSST patch
        pixel_scale=0.168, # HSC pixel scale
        mag_limit=mag_limit,
        mag_limit_band=mag_limit_band,
    )

    x0 = bbox.beginX
    y0 = bbox.beginY
    xcoords = ssp.x + x0
    ycoords = ssp.y + y0
    radec_coords = wcs.pixelToSkyArray(xcoords, ycoords, degrees=True)

    # Add extinction
    # see https://github.com/rubin-dp0/delegate-contributions-dp02/blob/db7d06ba6203faa15732c5368ff6e52ea53c5796/MWhalo_density/Milky_Way_halo_density.ipynb#L64

    # set the A_lamba/E(B-V) values for the six ugrizy LSST filters
    band_a_ebv = np.array([4.812, 3.643, 2.699, 2.063, 1.578, 1.313])

    coords = SkyCoord(radec_coords[0], radec_coords[1],
                      unit='deg', frame='icrs')

    sfd = SFDQuery()
    ebvvec = sfd(coords)

    u_ext = ssp.mags['LSST_u'] + ebvvec*band_a_ebv[0]
    g_ext = ssp.mags['LSST_g'] + ebvvec*band_a_ebv[1]
    r_ext = ssp.mags['LSST_r'] + ebvvec*band_a_ebv[2]
    i_ext = ssp.mags['LSST_i'] + ebvvec*band_a_ebv[3]
    z_ext = ssp.mags['LSST_z'] + ebvvec*band_a_ebv[4]

    if band_for_injection == 'g':
        mag_for_injection = g_ext
    elif band_for_injection == 'r':
        mag_for_injection = r_ext
    elif band_for_injection == 'i':
        mag_for_injection = i_ext
    else:
        mag_for_injection = i_ext

    dist_column = np.repeat(dist.value, len(ssp.mags))

    cat = Table({'injection_id': np.arange(len(ssp.mags)),
                 'ra': radec_coords[0], 'dec': radec_coords[1],
                 'source_type': ['DeltaFunction']*len(ssp.mags),
                 'distance': dist_column,
                 'g_mag': g_ext, 'r_mag': r_ext, 'i_mag': i_ext,
                 'mag': mag_for_injection,
                 })

    return cat


def massage_the_cat(cat_inp, injection_maglim, band_for_injection,
                    xcen, ycen, wcs, bbox,
                    r_scale=300.0, dist=2.0):
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
    xcen : `float`
        X-coordinate to center the dwarf on (between 0-4000)
    ycen : `float`
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

    # Calculate the total mag below the mag limit.
    # Reduce the original catalog to only stars above the magnitude limit.
    # Append the Sersic model for the remaining flux.

    band = band_for_injection+'_mag'
    mag_for_sersic = totmag_below_maglim(cat_inp[band], injection_maglim)

    # Replicate the magnitude column for the band you want to inject into:
    cat_inp.replace_column('mag', cat_inp[band])

    cat = cat_inp[cat_inp['mag'] <= injection_maglim]

    # Parameters for simulated dwarf:
    r_scale = r_scale*u.pc
    reff = (r_scale/1.68).to(u.kpc)
    pa = 0.0
    axis_ratio = 1.0
    # ra_sim = np.median(cat['ra'])
    # dec_sim = np.median(cat['dec'])

    # Append a single line for the "galaxy" that contains the unresolved flux:
    cat.add_row()

    x0 = bbox.beginX
    y0 = bbox.beginY
    xcoord = xcen + x0
    ycoord = ycen + y0
    radec_coords = wcs.pixelToSkyArray(xcoord, ycoord, degrees=True)

    cat[-1]['ra'] = radec_coords[0][0]
    cat[-1]['dec'] = radec_coords[1][0]
    cat[-1]['mag'] = mag_for_sersic
    cat[-1]['source_type'] = 'Sersic'

    # I think instead of including distance, I need to convert the radius
    # into a sky radius instead of physical. CHECK THIS!
    cat[-1]['distance'] = dist

    # We need Sersic index, position angle, axis ratio, and semimajor axis for the galaxy model,
    #   so create columns for these (with all stars set to some default values)
    semimajor_all = 0.0*cat['mag']
    semimajor_all[-1] = reff.value
    sersic_n_all = 0.0*cat['mag']
    sersic_n_all[-1] = 1.0
    pa_all = 0.0*cat['mag']
    axis_ratio_all = np.ones(len(cat['mag']))

    cat.add_columns([semimajor_all, sersic_n_all, pa_all, axis_ratio_all],
                    names=['r_eff', 'n', 'pa', 'axis_ratio'])

    return cat
