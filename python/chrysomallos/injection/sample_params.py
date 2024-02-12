# class GenerateParams():
#     distance_sampler:
import os

import fitsio
import numpy as np
import pandas as pd

from chrysomallos.utils import (
    get_coadd_dict,
    logger,
    mstar_from_absmag,
    rh_mv_to_sb,
    sb_mv_to_rh,
    sb_rh_to_mv,
)

__all__ = [
    "DwarfParamSampler",
]

OUT_ORDER = [
    "id",
    "tract",
    "patch",
    "x_cen",
    "y_cen",
    "ra",
    "dec",
    "distance",
    "m_v",
    "surface_brightness",
    "ellipticity",
    "theta",
    "age",
    "feh",
    "stellar_mass",
    "r_scale",
    "n",
    "random_seed_injection",
]


class DwarfParamSampler:
    """
    Class for generating dwarf parameters:
    """

    def __init__(self, config, coadd_dict=None):
        self.config = config
        self.coadd_dict = get_coadd_dict(coadd_dict=coadd_dict, config=self.config)
        # To Do: do we reading the bbox and adjust x/y_cen ranges?

    def run(self, write=True):
        sampling_config = self.config["sampling"]
        if sampling_config["random_seed_sampling"] is not None:
            np.random.seed(sampling_config["random_seed_sampling"])

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
        self.dwarf_param_frame["id"] = self.dwarf_param_frame.index
        # get mstar_from_absmag populate mstar
        self.dwarf_param_frame["stellar_mass"] = mstar_from_absmag(
            self.dwarf_param_frame["m_v"]
        )
        # get ra/dec
        ra, dec = self._populate_ra_dec()

        self.dwarf_param_frame["ra"] = ra
        self.dwarf_param_frame["dec"] = dec

        self.dwarf_param_frame = self.dwarf_param_frame[OUT_ORDER]
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
            rec_arr = self.dwarf_param_frame.to_records(index=False)
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

    def _populate_ra_dec(self):
        """
        Populate ra/dec from x/y_cen
        """
        # grab the wcs from the coadd
        first_band = self.config["pipelines"]["bands"][0]
        wcs = self.coadd_dict[first_band]["wcs"]

        ra, dec = wcs.pixelToSkyArray(
            self.dwarf_param_frame["x_cen"].values,
            self.dwarf_param_frame["y_cen"].values,
            degrees=True,
        )
        return ra, dec
