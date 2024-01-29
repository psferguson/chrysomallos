import unittest

import numpy as np
import scipy.stats as stats

from chrysomallos.synthpop import inverse_transform_sample, sersic_xy


class SpatialSamplerTest(unittest.TestCase):
    def test_inverse_transform_sample(self):
        """"""
        np.random.seed(123)
        samples = np.array(
            [
                0.50505051,
                -0.70707071,
                -0.90909091,
                0.1010101,
                0.50505051,
                -0.3030303,
                1.91919192,
                0.3030303,
                -0.1010101,
                -0.3030303,
            ]
        )
        vals = np.linspace(-10, 10, 100)
        pdf = stats.norm.pdf(vals, loc=0, scale=1)

        samples_out = np.round(inverse_transform_sample(vals, pdf, 10), 8)
        self.assertEqual(list(samples_out), list(samples))

    def test_sersic_xy(self):
        np.random.seed(123)
        samples = np.array(
            [
                [1835.31103058, 2049.05395149],
                [1982.99872288, 1974.12255281],
                [1897.71765116, 2007.9530405],
                [2211.28065438, 2017.37382499],
                [1750.27719196, 2036.3424305],
                [1986.75628041, 1964.57652989],
                [2296.92199151, 2132.74532411],
                [2131.27693085, 2052.42328724],
                [1805.06154694, 1991.48642921],
                [1836.96493766, 1992.82083895],
            ]
        )
        xy_kw = dict(
            num_stars=10,
            r_eff=1,
            n=1,
            theta=0.2,
            ellip=0.8,
            distance=2,
            xy_dim=4001,
            num_r_eff=10,
            dx=2000,
            dy=2000,
            pixel_scale=0.5,
            random_state=123,
        )
        samples_out = np.round(sersic_xy(**xy_kw).data, 8)
        self.assertEqual(samples_out.tolist(), samples.tolist())


if __name__ == "__main__":
    unittest.main()
