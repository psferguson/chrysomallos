# class GenerateParams():
#     distance_sampler:
import os

import fitsio
import numpy as np
import pandas as pd

from chrysomallos.utils import logger, rh_mv_to_sb, sb_mv_to_rh, sb_rh_to_mv

__all__ = [
    "DwarfParamSampler",
]


class DwarfParamSampler:
    """
    Class for generating dwarf parameters:
    """

    def __init__(self, config, seed=None):
        self.config = config
        self.seed = seed
        # To Do: do we reading the bbox and adjust x/y_cen ranges?

    def run(self, write=True):
        sampling_config = self.config["sampling"]
        if sampling_config["random_seed"] is not None:
            np.random.seed(sampling_config["random_seed"])

        if sampling_config["type"] == "grid":
            raise (NotImplementedError)

        elif sampling_config["type"] == "sample":
            n_dwarfs = sampling_config["n_dwarfs"]
            logger.info(f"sampling {n_dwarfs} dwarfs")
            output_dict = {key: [] for key in sampling_config["params"]}
            for param in output_dict.keys():
                param_val = sampling_config["params"][param]
                if (param in ["m_v", "surface_brightness", "r_scale"]) and (
                    param_val is None
                ):
                    output_dict[param] = np.ones(n_dwarfs) * np.nan
                    calc_param = param
                elif (param in ["ra", "dec", "stellar_mass"]) & (param_val is not None):
                    logger.info(f"sampling for param: {param} not yet implemented")
                    output_dict[param] = np.ones(n_dwarfs) * np.nan

                else:
                    output_dict[param] = self.sample(param_val, n_dwarfs)

        self.dwarf_param_frame = pd.DataFrame(output_dict)
        self.dwarf_param_frame[calc_param] = self._fill_mv_sb_r_scale(
            calc_param, self.dwarf_param_frame
        )
        self.dwarf_param_frame["tract"] = self.config["pipelines"]["tract"]
        self.dwarf_param_frame["patch"] = self.config["pipelines"]["patch"]
        self.dwarf_param_frame["dwarf_id"] = self.dwarf_param_frame.index
        # To Do:  mstar_from_absmag populate mstar/ra/dec
        if write:
            logger.info("saving generated params")
            self.write_param_file()
        return self.dwarf_param_frame

    def sample(self, param, n_dwarfs):
        """
        If passed a number return array of that value length n_dwarfs
        If passed list with 3 elements and last element is linear
        """
        if isinstance(param, float) | isinstance(param, int):
            data = np.ones(n_dwarfs) * param
        elif isinstance(param, list):
            if param[-1] == "linear":
                data = self.linear(param[0], param[1], n_dwarfs)
        elif param is None:
            data = np.ones(n_dwarfs) * np.nan
        else:
            raise Exception(f"Sample failed for {param}")

        return data

    def linear(self, min_val, max_val, size):
        if max_val <= min_val:
            raise Exception(f"Max ({max_val}) must be greater than Min ({min_val})")
        return np.random.uniform(min_val, max_val, size)

    def write_param_file(self):
        filename = self.config["sampling"]["output_file"]
        ext = os.path.splitext(filename)[1]
        if ext == ".csv":
            self.dwarf_param_frame.to_csv(filename, index=False)
        elif ext == ".fits":
            rec_arr = self.dwarf_param_frame.to_records()
            fitsio.write(filename, rec_arr, clobber=True)
        else:
            raise Exception(f"bad filetype {ext}")
        return 2

    def _fill_mv_sb_r_scale(self, calc_param, df):
        if calc_param == "m_v":
            vals = sb_rh_to_mv(
                sb=df["surface_brightness"],
                rh=df["surface_scale"],
                distance=df["distance"],
            )
        elif calc_param == "r_scale":
            vals = sb_mv_to_rh(
                sb=df["surface_brightness"], M_v=df["m_v"], distance=df["distance"]
            )
        elif calc_param == "surface_brightness":
            vals = rh_mv_to_sb(
                M_v=df["m_v"], rh=df["surface_scale"], distance=df["distance"]
            )
        else:
            raise Exception(f"unknown calc param: {calc_param}")
        return vals
