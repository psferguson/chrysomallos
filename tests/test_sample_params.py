import unittest

import lsst.afw.geom as afwGeom
import lsst.geom as geom
import numpy as np
from numpy.testing import assert_array_almost_equal

from chrysomallos import DwarfParamSampler


class TestDwarfParamSampler(unittest.TestCase):
    def setUp(self):
        self.r_scale_vals = np.array(
            [
                77.25481124,
                193.21962584,
                85.41326808,
                59.36009232,
                92.21115948,
                34.73158744,
                51.32471657,
                65.81732571,
                83.73329397,
                46.56359376,
            ]
        )

        self.config = {
            "pipelines": {
                "bands": ["g", "r", "i"],
                "input_collections": [
                    "HSC/runs/RC2/w_2023_32/DM-40356/20230819T003257Z"
                ],
                "patch": 3,
                "repo": "~/myprojects/rubin/dwarfFinder/deepCoadd_repo",
                "skymap": "hsc_rings_v1",
                "tract": 9615,
            },
            "sampling": {
                "n_dwarfs": 10,
                "params": {
                    "age": 10,
                    "dec": None,
                    "dist": [2, 3, "linear"],
                    "ellip": [0, 0.5, "linear"],
                    "feh": -2.0,
                    "m_v": [-9, -5, "linear"],
                    "r_scale": None,
                    "ra": None,
                    "n_sersic": 1,
                    "theta": [0, 180, "linear"],
                    "stellar_mass": None,
                    "surface_brightness": [24, 27, "linear"],
                    "x_cen": [300, 3800, "linear"],
                    "y_cen": [300, 3800, "linear"],
                    "random_seed_injection": None,
                },
                "output_file": "temp_params.csv",
                "output_directory": "temp_dir",
                "random_seed_sampling": 69,
                "type": "sample",
            },
        }
        self.coadd_dict = self.create_fake_coadd_dict()
        self.sampler = DwarfParamSampler(self.config, self.coadd_dict)

    def create_fake_coadd_dict(self):
        return {
            "g": {
                "image": [],
                "wcs": afwGeom.makeSkyWcs(
                    crpix=geom.Point2D(17999, 17999),
                    crval=geom.SpherePoint(110.0, 39.3, geom.degrees),
                    cdMatrix=afwGeom.makeCdMatrix(scale=0.175 * geom.arcseconds),
                ),
                "bbox": geom.Box2I(
                    corner=geom.Point2I(11900, 0),
                    dimensions=geom.Extent2I(4200, 4100),
                ),
                "psf": [],
                "photo_calib": [],
            }
        }

    def test_regular_sample(self):
        dwarf_df, _ = self.sampler.run()
        self.assertEqual(dwarf_df.shape, (10, 18))
        assert_array_almost_equal(dwarf_df["r_scale"].values, self.r_scale_vals)

    def test_bad_sample_param(self):
        old_config = self.config.copy()
        self.config["sampling"]["params"]["test"] = [0, 10, "linear"]
        sampler = DwarfParamSampler(self.config, self.coadd_dict)
        with self.assertRaises(AssertionError):
            sampler.run()
        self.config = old_config
