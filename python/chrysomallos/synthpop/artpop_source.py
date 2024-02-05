# Third-party
import numpy as np
from artpop import MIST_PATH, Source
from artpop.stars import MISTSSP
from artpop.util import check_units, check_xy_dim
from astropy import units as u
from scipy.special import gamma, gammaincinv

from chrysomallos.synthpop.spatial_sampler import sersic_xy

__all__ = [
    "SersicSPChrysomallos",
    "MISTSersicSSPChrysomallos",
]


class SersicSPChrysomallos(Source):
    """
    Stellar population with a Sersic spatial distribution.

    Parameters
    ----------
    sp : `~artpop.stars.populations.StellarPopulation`
        A stellar population object.
    r_eff : float or `~astropy.units.Quantity`
        Effective radius of the source. If a float is given, the units are
        assumed to be `kpc`. Must be greater than zero.
    n : float
        Sersic index. Must be greater than zero.
    theta : float or `~astropy.units.Quantity`
        Rotation angle, counterclockwise from the positive x-axis. If a float
        is given, the units are assumed to be `degree`.
    ellip : float
        Ellipticity.
    xy_dim : int or list-like
        Dimensions of the mock image in xy coordinates. If int is given,
        will make the x and y dimensions the same.
    pixel_scale : float or `~astropy.units.Quantity`
        The pixel scale of the mock image. If a float is given, the units will
        be assumed to be `arcsec / pixels`. Default is `0.2 arcsec / pixel`.
    num_r_eff : float, optional
        Number of r_eff to sample positions within. This parameter is needed
        because the current Sersic sampling function samples from within a
        discrete grid. Default is 10.
    dx : float, optional
        Shift from center of image in the x direction.
    dy : float, optional
        Shift from center of image in the y direction.
    labels : list-like, optional
        Labels for the stars. For example, EEP values (int or float) or name
        of evolutionary phase (str).
    """

    def __init__(
        self,
        sp,
        r_eff,
        n,
        theta,
        ellip,
        xy_dim,
        pixel_scale,
        num_r_eff=10,
        dx=0,
        dy=0,
        labels=None,
    ):
        self.sp = sp
        self.mag_limit = sp.mag_limit
        self.mag_limit_band = sp.mag_limit_band
        self.smooth_model = None

        if self.mag_limit is not None and sp.frac_num_sampled < 1.0:
            _r_eff = check_units(r_eff, "kpc").to("Mpc").value
            _theta = check_units(theta, "deg").to("radian").value
            _distance = check_units(sp.distance, "Mpc").to("Mpc").value
            _pixel_scale = check_units(pixel_scale, u.arcsec / u.pixel)

            if _r_eff <= 0:
                raise Exception("Effective radius must be greater than zero.")

            xy_dim = check_xy_dim(xy_dim)
            x_0, y_0 = xy_dim // 2
            x_0 += dx
            y_0 += dy
            self.n = n
            self.ellip = ellip
            self.r_sky = np.arctan2(_r_eff, _distance) * u.radian.to("arcsec")
            self.r_sky *= u.arcsec
            r_pix = self.r_sky.to("pixel", u.pixel_scale(_pixel_scale)).value

        self.xy_kw = dict(
            num_stars=sp.num_stars,
            r_eff=r_eff,
            n=n,
            theta=theta,
            ellip=ellip,
            distance=sp.distance,
            xy_dim=xy_dim,
            num_r_eff=num_r_eff,
            dx=x_0 + dx,
            dy=y_0 + dy,
            pixel_scale=pixel_scale,
            random_state=sp.rng,
        )

        _xy = sersic_xy(**self.xy_kw)

        super(SersicSPChrysomallos, self).__init__(
            _xy, sp.mag_table, xy_dim, pixel_scale, labels
        )

    def mag_to_image_amplitude(self, m_tot, zpt):
        """
        Convert total magnitude into amplitude parameter for the smooth model.

        Parameters
        ----------
        m_tot : float
            Total magnitude in the smooth component of the system.
        zpt : float
            Photometric zero point.

        Returns
        -------
        mu_e : float
            Surface brightness at the effective radius of the Sersic
            distribution in mags per square arcsec.
        amplitude : float
            Amplitude parameter for the smooth model in image flux units.
        param_name : str
            Name of amplitude parameter (needed to set its value when
            generating the smooth model).
        """
        param_name = "amplitude"
        b_n = gammaincinv(2.0 * self.n, 0.5)
        f_n = gamma(2 * self.n) * self.n * np.exp(b_n) / b_n ** (2 * self.n)
        r_circ = self.r_sky * np.sqrt(1 - self.ellip)
        area = np.pi * r_circ.to("arcsec").value ** 2
        mu_e = m_tot + 2.5 * np.log10(2 * area) + 2.5 * np.log10(f_n)
        amplitude = 10 ** (0.4 * (zpt - mu_e)) * self.pixel_scale.value**2
        return mu_e, amplitude, param_name


def _check_label_type(ssp, label_type, has_phases=True):
    if label_type is not None:
        if label_type == "phases" and has_phases:
            labels = ssp.get_star_phases()
        elif hasattr(ssp, label_type):
            labels = getattr(ssp, label_type)
        else:
            raise Exception(f"{label_type} is not a valid label type.")
    else:
        labels = None
    return labels


class MISTSersicSSPChrysomallos(SersicSPChrysomallos):
    """
    MIST simple stellar population with a Sersic spatial distribution. This
    is a convenience class that combines `~artpop.space.sersic_xy` and
    `~artpop.stars.MISTSSP` to make a `~artpop.source.Source` object.

    .. note::
        You must give `total_mass` *or* `num_stars`.

    Parameters
    ----------
    log_age : float
        Log (base 10) of the simple stellar population age in years.
    feh : float
        Metallicity [Fe/H] of the simple stellar population.
    phot_system : str or list-like
        Name of the photometric system(s).
    r_eff : float or `~astropy.units.Quantity`
        Effective radius of the source. If a float is given, the units are
        assumed to be `kpc`. Must be greater than zero.
    n : float
        Sersic index. Must be greater than zero.
    theta : float or `~astropy.units.Quantity`
        Rotation angle, counterclockwise from the positive x-axis. If a float
        is given, the units are assumed to be `degree`.
    ellip : float
        Ellipticity defined as `1 - b/a`, where `b` is the semi-minor axis
        and `a` is the semi-major axis.
    distance : float or `~astropy.units.Quantity`
        Distance to source. If float is given, the units are assumed
        to be `Mpc`.
    xy_dim : int or list-like
        Dimensions of the mock image in xy coordinates. If int is given,
        will make the x and y dimensions the same.
    pixel_scale : float or `~astropy.units.Quantity`
        The pixel scale of the mock image. If a float is given, the units will
        be assumed to be `arcsec / pixels`. Default is `0.2 arcsec / pixel`.
    num_stars : int or `None`
        Number of stars in source. If `None`, then must give `total_mass`.
    total_mass : float or `None`
        Stellar mass of the source. If `None`, then must give `num_stars`. This
        mass accounts for stellar remnants, so the actual sampled mass will be
        less than the given value.
    a_lam : float or dict, optional
        Magnitude(s) of extinction. If float, the same extinction will be applied
        to all bands. If dict, the keys must be the same as the observational filters.
    add_remnants : bool, optional
        If True (default), apply scaling factor to total mass to account for
        stellar remnants in the form of white dwarfs, neutron stars,
        and black holes.
    mag_limit : float, optional
        Only sample individual stars that are brighter than this magnitude. All
        fainter stars will be combined into an integrated component. Otherwise,
        all stars in the population will be sampled. You must also give the
        `mag_limit_band` if you use this parameter.
    mag_limit_band : str, optional
        Bandpass of the limiting magnitude. You must give this parameter if
        you use the `mag_limit` parameter.
    imf : str, optional
        The initial stellar mass function. Default is `'kroupa'`.
    imf_kw : dict, optional
        Optional keyword arguments for sampling the stellar mass function.
    mist_path : str, optional
        Path to MIST isochrone grids. Use this if you want to use a different
        path from the `MIST_PATH` environment variable.
    num_r_eff : float, optional
        Number of r_eff to sample positions within. This parameter is needed
        because the current Sersic sampling function samples from within a
        discrete grid. Default is 10.
    mass_tolerance : float, optional
        Tolerance in the fractional difference between the input mass and the
        final mass of the population. The parameter is only used when
        `total_mass` is given.
    dx : float, optional
        Shift from center of image in the x direction.
    dy : float, optional
        Shift from center of image in the y direction.
    label_type : str, optional
        If not None, the type of labels to pass to the `~artpop.source.Source`
        object. Can be `phases` or an attribute of `~artpop.stars.MISTSSP`.
    random_state : `None`, int, list of ints, or `~numpy.random.RandomState`
        If `None`, return the `~numpy.random.RandomState` singleton used by
        ``numpy.random``. If `int`, return a new `~numpy.random.RandomState`
        instance seeded with the `int`.  If `~numpy.random.RandomState`,
        return it. Otherwise raise ``ValueError``.
    """

    def __init__(
        self,
        log_age,
        feh,
        phot_system,
        r_eff,
        n,
        theta,
        ellip,
        distance,
        xy_dim,
        pixel_scale,
        num_stars=None,
        total_mass=None,
        a_lam=0.0,
        add_remnants=True,
        mag_limit=None,
        mag_limit_band=None,
        imf="kroupa",
        imf_kw=None,
        mist_path=MIST_PATH,
        num_r_eff=10,
        mass_tolerance=0.01,
        dx=0,
        dy=0,
        label_type=None,
        random_state=None,
    ):
        self.ssp_kw = dict(
            log_age=log_age,
            feh=feh,
            phot_system=phot_system,
            distance=distance,
            a_lam=a_lam,
            total_mass=total_mass,
            num_stars=num_stars,
            imf=imf,
            mist_path=mist_path,
            imf_kw=imf_kw,
            random_state=random_state,
            mag_limit=mag_limit,
            mag_limit_band=mag_limit_band,
            add_remnants=add_remnants,
            mass_tolerance=mass_tolerance,
        )

        ssp = MISTSSP(**self.ssp_kw)
        labels = _check_label_type(ssp, label_type)

        super(MISTSersicSSPChrysomallos, self).__init__(
            sp=ssp,
            r_eff=r_eff,
            n=n,
            theta=theta,
            ellip=ellip,
            xy_dim=xy_dim,
            pixel_scale=pixel_scale,
            num_r_eff=num_r_eff,
            dx=dx,
            dy=dy,
            labels=labels,
        )
