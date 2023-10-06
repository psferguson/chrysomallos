import numpy as np
import scipy.stats
from astropy.modeling.models import Sersic1D



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
        points=self.sample_position(radii, self.e, self.theta)
        return points
        
class SersicSampler(SpatialSampler):
    """ Sample from interpolated function. """
    def __init__(self, r_eff,n,theta,e, **kwargs):
        """ Interpolation sampler.

        Parameters
        ----------
        xvals : x-values of interpolation
        yvals : y-values of interpolation
        kwargs : passed to baseclass

        Returns
        -------
        sampler
        """
        self.r_eff = r_eff
        self.n = n
        self.theta = theta
        self.e = e
        self.interp = Sersic1D(amplitude=1, r_eff = r_eff, n=n)
        self.interp.xmin=0
        self.interp.xmax= self.r_eff * 10

    def sample_radii(self, size, nsteps=1e5):
        xvals = np.linspace(self.interp.xmin, self.interp.xmax, int(nsteps))
        pdf = self.interp(xvals) * 2 * np.pi * xvals 
        return inverse_transform_sample(xvals, pdf, size=size)
    
    def sample_position(self, radii, e, theta):
        angle = 2 * np.pi * np.random.uniform(0, 1, size=len(radii))
        x = radii * np.cos(angle)
        y = radii * np.sin(angle)
        y *= (1 - e) # make an ellipse
        
        # rotate points by theta
        rot_mat=np.array([[np.cos(theta), -np.sin(theta)],
                         [np.sin(theta), np.cos(theta)]])
        x_rot,y_rot = rot_mat.dot(np.vstack([x,y]))
        return x_rot,y_rot
        