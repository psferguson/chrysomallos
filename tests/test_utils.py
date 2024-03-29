import unittest

import lsst.afw.geom as afwGeom
import lsst.geom as geom
import numpy as np

from chrysomallos.injection import adopt_a_cat
from chrysomallos.utils import dist_to_dmod, sdss_g_to_V, totmag


class MVTest(unittest.TestCase):
    def test_mv_application(self):
        mv_vals = np.arange(-10, -8, 0.5)

        # Confirm that passing each of these M_V values to the
        # sythetic catalog generation code returns the correct
        # stellar mass.

        # adopt_a_cat needs a WCS and BBox
        bbox_fake = geom.Box2I(
            corner=geom.Point2I(11900, 0),
            dimensions=geom.Extent2I(4200, 4100),
        )
        wcs_fake = afwGeom.makeSkyWcs(
            crpix=geom.Point2D(17999, 17999),
            crval=geom.SpherePoint(110.0, 39.3, geom.degrees),
            cdMatrix=afwGeom.makeCdMatrix(scale=0.175 * geom.arcseconds),
        )

        dist = 1.0  # Distance in Mpc
        MV_out = []
        for m_v in mv_vals:
            cat = adopt_a_cat(
                wcs_fake,
                bbox_fake,
                age=10.0,
                feh=-2.0,
                dist=1.0,
                r_scale=100.0,
                ellip=0,
                theta=0,
                n=1,
                m_v=m_v,
                mag_limit=36.0,
                mag_limit_band="LSST_g",
                random_seed=None,
            )
            vmags = sdss_g_to_V(cat["g_mag"], cat["r_mag"])
            MV_out = totmag(vmags) - dist_to_dmod(dist * 1.0e6)
            self.assertAlmostEqual(m_v, MV_out, delta=0.2)


if __name__ == "__main__":
    unittest.main()
