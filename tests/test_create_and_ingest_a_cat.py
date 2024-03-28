import unittest

# import lsst.afw.geom as afwGeom
# import lsst.geom as geom
# import numpy as np
# from numpy.testing import assert_array_almost_equal
from test_sample_params import create_fake_coadd_dict

from chrysomallos import CreateDwarfInjectionCatalog


class TestCreateDwarfInjectionCatalog(unittest.TestCase):
    def setUp(self):
        self.config = {}
        self.dwarf_params_frame = {}
        self.coadd_dict = create_fake_coadd_dict()
        self.creator = CreateDwarfInjectionCatalog()

    def test_regular_sample(self):
        import pdb

        pdb.set_trace()
