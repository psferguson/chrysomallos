import unittest

import lsst.afw.geom as afwGeom
import lsst.geom as geom
from astropy.table import Table

from chrysomallos.injection.adopt_a_cat import adopt_a_cat, massage_the_cat


class TestAdoptACat(unittest.TestCase):
    def setUp(self):
        self.bbox_fake = geom.Box2I(
            corner=geom.Point2I(11900, 0),
            dimensions=geom.Extent2I(4200, 4100),
        )
        self.wcs_fake = afwGeom.makeSkyWcs(
            crpix=geom.Point2D(17999, 17999),
            crval=geom.SpherePoint(110.0, 39.3, geom.degrees),
            cdMatrix=afwGeom.makeCdMatrix(scale=0.175 * geom.arcseconds),
        )

    def test_adopt_a_cat(self):
        cat = adopt_a_cat(self.wcs_fake, self.bbox_fake)
        self.assertIsInstance(cat, Table)
        expected_columns = [
            "injection_id",
            "ra",
            "dec",
            "source_type",
            "dist",
            "g_mag",
        ]
        self.assertTrue(all(column in cat.columns for column in expected_columns))


class TestMassageTheCat(unittest.TestCase):
    def setUp(self):
        self.cat_inp = Table(
            {"ra": [10.0, 20.0], 
             "dec": [30.0, 40.0], 
             "LSST_g_mag": [15.0, 16.0],
             "source_type": ["DeltaFunction", "DeltaFunction"],
             "dist": [2.0, 2.0],
            }
        )
        self.mag_limit = 20.0
        self.band_for_injection = "LSST_g"
        self.bbox = geom.Box2I(
            corner=geom.Point2I(11900, 0),
            dimensions=geom.Extent2I(4200, 4100),
        )  
        self.wcs = afwGeom.makeSkyWcs(
            crpix=geom.Point2D(17999, 17999),
            crval=geom.SpherePoint(110.0, 39.3, geom.degrees),
            cdMatrix=afwGeom.makeCdMatrix(scale=0.175 * geom.arcseconds),
        )  
        self.x_cen = 15.0
        self.y_cen = 35.0

    def test_massage_the_cat(self):
        cat = massage_the_cat(
            self.cat_inp,
            self.mag_limit,
            self.band_for_injection,
            self.wcs,
            self.bbox,
            self.x_cen,
            self.y_cen,
        )
        self.assertIsInstance(cat, Table)
        expected_columns = ["ra", "dec", "mag"]
        self.assertTrue(all(column in cat.columns for column in expected_columns))


if __name__ == "__main__":
    unittest.main()
