# class GenerateParams():
#     distance_sampler:
import numpy as np
import pandas as pd

from chrysomallos.utils import mstar_from_absmag, sb_mv_to_rh


class DwarfParamSampler:
    """
    Class for generating dwarf parameters:
    """

    def __init__(self, config, seed=None):
        self.config = config
        self.seed = seed

        # if sampler_type == "random_uniform":
        #     self.sampler = self.random_uniform_sampler

    def run(self):
        sampling_config = self.config["sampling"]
        n_dwarfs = sampling_config["n_dwarfs"]
        # output_dict = [key:[] for key in sampling_config["params"]]
        # for param in

    def random_uniform_sampler(self, stamp_config):
        distances = np.random.uniform(
            stamp_config.dist_range[0],
            stamp_config.dist_range[1],
            stamp_config.n_dwarfs,
        )
        Mvs = np.random.uniform(
            stamp_config.Mv_range[0], stamp_config.Mv_range[1], stamp_config.n_dwarfs
        )
        sbs = np.random.uniform(
            stamp_config.sb_range[0], stamp_config.sb_range[1], stamp_config.n_dwarfs
        )
        r_scales = sb_mv_to_rh(sb=sbs, M_v=Mvs, distance=distances * 1e6)

        ellips = np.random.uniform(
            stamp_config.ellip_range[0],
            stamp_config.ellip_range[1],
            stamp_config.n_dwarfs,
        )
        # assuming all cens tract coordinates
        x_cens = np.random.uniform(
            stamp_config.x_cen_range[0],
            stamp_config.x_cen_range[1],
            stamp_config.n_dwarfs,
        )

        y_cens = np.random.uniform(
            stamp_config.y_cen_range[0],
            stamp_config.y_cen_range[1],
            stamp_config.n_dwarfs,
        )

        # tract, patch, id, dist, x_cen, y_cen, r_scale, sb, mag_limit, m_v, mass, random_seed
        rng = np.random.default_rng()
        seeds = rng.integers(-2147483647, high=2147483646, size=stamp_config.n_dwarfs)
        param_frame = pd.DataFrame(
            {
                "tract": stamp_config.tract,
                "patch": stamp_config.patch,
                "id": np.arange(stamp_config.n_dwarfs),
                "dist": distances,
                "x_cen": x_cens,
                "y_cen": y_cens,
                "r_scale": r_scales,
                "sb": sbs,
                "Mv": Mvs,
                "mass": mstar_from_absmag(Mvs),
                "mag_limit": stamp_config.mag_lim,
                "random_seed": seeds,
            }
        )
        return param_frame
