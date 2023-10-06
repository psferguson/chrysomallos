import numpy as np
import scipy.stats
from astropy.modeling.models import Sersic1D
from artpop.util import check_units, check_xy_dim
import astropy.units as u

def inverse_transform_sample(vals, pdf, size):
    """ Perform inverse transform sampling

    Parameters
    ----------
    vals: value at which pdf is measured
    pdf : pdf value
    size : number of stars to sample

    Returns
    -------
    samples : samples of vals
    """
    
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    fn = scipy.interpolate.interp1d(cdf, list(range(0, len(cdf))),
                                    bounds_error=False,
                                    fill_value = 0.0)
    x_new = scipy.stats.uniform.rvs(size=np.rint(size).astype(int))
    index = np.rint(fn(x_new)).astype(int)
    return vals[index]

class SpatialSampler(object):
    def __init__(self, *args, **kwargs):
        pass

    def pdf(self, values):
        return self._pdf(values)

    def sample_radii(self, size):
        pass
    def sample_position(self):
        pass 
    def sample(self, size):
        radii=self.sample_radii(size)
        points=self.sample_position(radii, self.ellip, self.theta)
        return points
        
class SersicSampler(SpatialSampler):
    """Sample from the Sersic profile.
        
    Attributes
    ----------
    r_eff : float
        The effective radius of the Sersic profile.
    n : float
        The Sersic index.
    ellip : float, optional
        The ellipticity of the profile.
    theta : float, optional
        The rotation angle of the profile in radians.
    x_0 : float, optional
        The x-coordinate of the profile's center.
    y_0 : float, optional
        The y-coordinate of the profile's center.
    amplitude : float, optional
        Amplitude of the Sersic profile.
    interp : Sersic1D object
        The Sersic profile function.
    """
    def __init__(self, r_eff, n, ellip=0, theta=0, x_0=0, y_0=0, amplitude=1,
                 num_r_eff=10, **kwargs):
        """ Initialize the SersicSampler.

        Returns
        -------
        sampler
        """
        self.r_eff = r_eff
        self.n = n
        self.theta = theta
        self.ellip = ellip
        self.interp = Sersic1D(amplitude=amplitude, 
                               r_eff = r_eff, 
                               n=n)
        self.interp.xmin=0
        self.interp.xmax= self.r_eff * num_r_eff
        self.x_0=x_0
        self.y_0=y_0

    def sample_radii(self, size, nsteps=1e5):
        """Sample radii from the Sersic profile.
        
        Parameters
        ----------
        size : int
            Number of samples to generate.
        nsteps : float, optional
            Number of steps for the linspace function (default is 1e5).
        
        Returns
        -------
        ndarray
            Array of sampled radii.
        """
        xvals = np.linspace(self.interp.xmin, self.interp.xmax, int(nsteps))
        pdf = self.interp(xvals) * 2 * np.pi * xvals 
        return inverse_transform_sample(xvals, pdf, size=size)
    
    def sample_position(self, radii, e, theta):
        """Sample positions (x, y) based on given radii, ellipticity, and rotation.
        
        Parameters
        ----------
        radii : ndarray
            Array of radii values.
        e : float
            The ellipticity of the profile.
        theta : float
            The rotation angle of the profile in radians.
        
        Returns
        -------
        tuple of ndarrays
            Tuple containing the x and y coordinates of the sampled positions.
        """
        angle = 2 * np.pi * np.random.uniform(0, 1, size=len(radii))
        x = radii * np.cos(angle)
        y = radii * np.sin(angle)
        y *= (1 - e) # make an ellipse
        
        # rotate points by theta
        rot_mat=np.array([[np.cos(theta), -np.sin(theta)],
                         [np.sin(theta), np.cos(theta)]])
        x_rot,y_rot = rot_mat.dot(np.vstack([x,y]))
        return x_rot + self.x_0, y_rot + self.y_0

def sersic_xy(num_stars, distance, xy_dim, pixel_scale, r_eff=500*u.pc, n=1.0,
              theta=45*u.deg, ellip=0.3, num_r_eff=10, dx=0, dy=0,
              drop_outside=False, random_state=None):
    """
    Sample xy positions from a two-dimensional Sersic distribution.


    .. note::
        Unlike artpop this uses an inverse transform sampler to be more memory
        efficient. The api is built to be the same as artpop.sersic_xy

    Parameters
    ----------
    num_stars : int
        Number of stars (i.e., positions) to sample.
    distance : float or `~astropy.units.Quantity`
        Distance to source. If float is given, the units are assumed
        to be `~astropy.units.Mpc`.
    xy_dim : list-like
        Dimensions of the mock image in xy coordinates. If int is given,
        will make the x and y dimensions the same.
    pixel_scale : float or `~astropy.units.Quantity`
        The pixel scale of the mock image. If a float is given, the units will
        be assumed to be `~astropy.units.arcsec` per `~astropy.units.pixels`.
    r_eff : float or `~astropy.units.Quantity`, optional
        Effective radius of the source. If a float is given, the units are
        assumed to be `~astropy.units.kpc`. Must be greater than zero.
    n : float, optional
        Sersic index. Must be greater than zero.
    theta : float or `~astropy.units.Quantity`, optional
        Rotation angle, counterclockwise from the positive x-axis. If a float
        is given, the units are assumed to be `degree`.
    ellip : float, optional
        Ellipticity defined as `1 - b/a`, where `b` is the semi-minor axis
        and `a` is the semi-major axis.
    num_r_eff : float, optional
        Number of r_eff to sample positions within. This parameter is needed
        because the current Sersic sampling function samples from within a
        discrete grid. Default is 10.
    dx : float, optional
        Shift from center of image in pixels in the x direction.
    dy : float, optional
        Shift from center of image in pixels in the y direction.
    drop_outside : bool, optional
        If True, drop all stars that fall outside of the image. In this case,
        the returned `~numpy.ndarray` will not be masked.

    Returns
    -------
    xy : `~numpy.ma.MaskedArray`
        Masked numpy array of xy positions. Positions that fall outside the
        mock image are masked.
    """
    if n <= 0:
        raise Exception('Sersic index n must be greater than zero.')

    xy_dim = check_xy_dim(xy_dim)       

    r_eff = check_units(r_eff, 'kpc').to('Mpc').value
    theta = check_units(theta, 'deg').to('radian').value
    distance = check_units(distance, 'Mpc').to('Mpc').value
    pixel_scale = check_units(pixel_scale, u.arcsec / u.pixel)

    if r_eff <= 0:
        raise Exception('Effective radius must be greater than zero.')

    r_pix = np.arctan2(r_eff, distance) * u.radian.to('arcsec') * u.arcsec
    r_pix = r_pix.to('pixel', u.pixel_scale(pixel_scale)).value
    sample_dim = 2 * np.ceil(r_pix * num_r_eff).astype(int) + 1
    
    sampler = SersicSampler(x_0=dx, y_0=dy, amplitude=1, r_eff=r_pix,
                            n=n, ellip=ellip, theta=theta, num_r_eff=num_r_eff)
    
    x,y = sampler.sample(num_stars)
    
    xy = np.vstack([x, y]).T
    xy = np.ma.masked_array(xy, mask=np.zeros_like(xy, dtype=bool))

    if xy_dim[0] < sample_dim or xy_dim[1] < sample_dim:
        outside_image = x < 0
        outside_image |= x > xy_dim[0] - 1
        outside_image |= y < 0
        outside_image |= y > xy_dim[1] - 1

        if outside_image.sum() > 0:
            msg = '{} stars outside the image'.format(outside_image.sum())
            print(msg)
            xy.mask = np.column_stack((outside_image, outside_image))

            if drop_outside:
                xy = xy[~outside_image].data
    return xy