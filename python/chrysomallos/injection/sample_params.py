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
    "dist",
    "m_v",
    "surface_brightness",
    "ellip",
    "theta",
    "age",
    "feh",
    "stellar_mass",
    "r_scale",
    "n_sersic",
    "n",
    "random_seed_injection",
    "dwarf_generation_id",
]

PARAMS_NOT_SAMPLED = [
    "ra",
    "dec",
    "stellar_mass",
    "random_seed_injection",
]


class DwarfParamSampler:
    """
    Samples parameters for dwarf galaxies based on configuration settings.
    This includes sampling galaxy properties and calculating dependent parameters.
    """

    def __init__(self, config, coadd_dict=None):
        """
        Initializes the sampler with a configuration dictionary and optional coaddition data.

        Parameters:
        - config: Configuration dictionary with sampling and pipeline settings.
        - coadd_dict: Optional dictionary containing coadded image data.
        """
        self.config = config
        self.coadd_dict = get_coadd_dict(coadd_dict=coadd_dict, config=self.config)
        os.makedirs(self.config["sampling"]["output_directory"], exist_ok=True)
        # TODO: do we want to read the bbox and adjust x/y_cen ranges?

    def run(self, write=True):
        """
        Executes the sampling process and optionally writes output to a file.

        Parameters:
        - write: Boolean indicating whether to write the sampled parameters to a file.

        Returns:
        - DataFrame containing sampled parameters for dwarf galaxies.
        """

        sampling_config = self.config["sampling"]
        grid_bool = False
        if sampling_config["random_seed_sampling"] is not None:
            np.random.seed(sampling_config["random_seed_sampling"])

        if sampling_config["type"] == "grid":
            raise (NotImplementedError)

        elif sampling_config["type"] == "sample":
            n_dwarfs = sampling_config["n_dwarfs"]
            # logger.info(f"sampling {n_dwarfs} dwarfs")
            output_dict = {key: [] for key in sampling_config["params"]}
            for param in output_dict.keys():
                param_val = sampling_config["params"][param]
                if (param in ["m_v", "surface_brightness", "r_scale"]) and (
                    param_val is None
                ):
                    output_dict[param] = np.ones(n_dwarfs) * np.nan
                    calc_param = param
                elif (param in PARAMS_NOT_SAMPLED) & (param_val is not None):
                    logger.info(f"sampling for param: {param} not yet implemented")
                    output_dict[param] = np.ones(n_dwarfs) * np.nan
                    # for any other parameter raise exception if is None or nan
                elif (param_val is None) & (param not in PARAMS_NOT_SAMPLED):
                    raise Exception(
                        f"param {param} is None or nan and cannot be sampled"
                    )
                elif param in ["x_cen", "y_cen"]:
                    if param_val[-1] == "grid":
                        grid_bool = True
                    else:
                        output_dict[param] = self.sample(param_val, n_dwarfs)
                else:
                    output_dict[param] = self.sample(param_val, n_dwarfs)
        if grid_bool:
            param_x_cen = sampling_config["params"]["x_cen"]
            param_y_cen = sampling_config["params"]["y_cen"]
            assert (param_x_cen[-1] == "grid") & (
                param_y_cen[-1] == "grid"
            ), "x/y_cen must be grid"
            param_val = param_x_cen[:-1] + param_y_cen[:-1] + ["grid"]
            grid_points = self.sample(param_val, n_dwarfs)
            output_dict["x_cen"] = grid_points[:, 0]
            output_dict["y_cen"] = grid_points[:, 1]

        self.dwarf_param_frame = pd.DataFrame(output_dict)
        self.dwarf_param_frame[calc_param] = self._fill_mv_sb_r_scale(
            calc_param, self.dwarf_param_frame
        )
        self.dwarf_param_frame["tract"] = self.config["pipelines"]["tract"]
        self.dwarf_param_frame["patch"] = self.config["pipelines"]["patch"]
        self.dwarf_param_frame["id"] = self.dwarf_param_frame.index
        self.dwarf_param_frame["dwarf_generation_id"] = self.config["sampling"][
            "generation_id"
        ]
        # get mstar_from_absmag populate mstar
        self.dwarf_param_frame["stellar_mass"] = mstar_from_absmag(
            self.dwarf_param_frame["m_v"]
        )
        # get ra/dec
        ra, dec = self._populate_ra_dec()

        self.dwarf_param_frame["ra"] = ra
        self.dwarf_param_frame["dec"] = dec
        self.dwarf_param_frame["n"] = self.dwarf_param_frame["n_sersic"]

        assert set(self.dwarf_param_frame.columns) == set(OUT_ORDER), (
            f"DataFrame columns mismatch. Expected: {set(OUT_ORDER)}"
            f", Found: {set(self.dwarf_param_frame.columns)}"
        )
        self.dwarf_param_frame = self.dwarf_param_frame[OUT_ORDER]
        if write:
            # logger.info("saving generated params")
            self.write_param_file()
        return self.dwarf_param_frame, self.coadd_dict

    def sample(self, param, n_dwarfs):
        """
        Samples a parameter based on its configuration.

        If passed a number return array of that value length n_dwarfs
        If passed list with 3 elements and last element is linear will
        uniformly sample between the first two elements.

        Parameters:
        - param: The parameter configuration to sample from.
        - n_dwarfs: The number of dwarf galaxies to sample parameters for.

        Returns:
        - An array of sampled values.
        """
        if isinstance(param, (float, int)):
            data = np.ones(n_dwarfs) * param
        elif isinstance(param, list):
            if param[-1] == "linear":
                data = self.linear(param[0], param[1], n_dwarfs)
            if param[-1] == "loglinear":
                data = self.loglinear(param[0], param[1], n_dwarfs)
            elif param[-1] == "grid":
                data = self.grid(param[0], param[1], param[2], param[3], n_dwarfs)
            elif all(isinstance(item, (float, int)) for item in param):
                data = np.zeros(n_dwarfs)
                for i in range(n_dwarfs):
                    data[i] = param[i % len(param)]
        elif param is None:
            data = np.ones(n_dwarfs) * np.nan
        else:
            raise Exception(f"Sample failed for {param}")

        return data

    def linear(self, min_val, max_val, size):
        """
        Generates linearly spaced samples between a minimum and maximum value.

        Parameters:
        - min_val: The minimum value in the range.
        - max_val: The maximum value in the range.
        - size: The number of samples to generate.

        Returns:
        - An array of linearly spaced samples.
        """
        if max_val <= min_val:
            raise Exception(f"Max ({max_val}) must be greater than Min ({min_val})")
        return np.random.uniform(min_val, max_val, size)

    def loglinear(self, min_val, max_val, size):
        """
        Generates linearly spaced samples between a minimum and maximum value.

        Parameters:
        - min_val: The minimum value in the range.
        - max_val: The maximum value in the range.
        - size: The number of samples to generate.

        Returns:
        - An array of linearly spaced samples.
        """
        if max_val <= min_val:
            raise Exception(f"Max ({max_val}) must be greater than Min ({min_val})")
        return 10 ** np.random.uniform(min_val, max_val, size)

    def grid(self, min_x, max_x, min_y, max_y, size):
        """
        Given number of points, generate a linear grid spanning x and y ranges.

        Parameters:
        min_x (float): Minimum x-coordinate.
        max_x (float): Maximum x-coordinate.
        min_y (float): Minimum y-coordinate.
        max_y (float): Maximum y-coordinate.
        size (int): Total number of points.

        Returns:
        np.ndarray: Array of shape (size, 2) containing the grid points.
        """
        assert size > 0, "size must be greater than 0"
        assert min_x < max_x, "min_x must be less than max_x"
        assert min_y < max_y, "min_y must be less than max_y"

        # Calculate the number of rows and columns
        n_rows = int(np.sqrt(size))
        n_cols = size // n_rows
        if n_rows * n_cols < size:
            n_cols += 1

        # Generate linearly spaced values within the given range for x and y
        x = np.linspace(min_x, max_x, n_cols)
        y = np.linspace(min_y, max_y, n_rows)

        # Create a meshgrid from these values
        xx, yy = np.meshgrid(x, y)

        # Combine the x and y coordinates into a single array of shape (n_rows*n_cols, 2)
        grid_points = np.c_[xx.ravel(), yy.ravel()]

        # If there are more points than needed, trim the excess
        if len(grid_points) > size:
            grid_points = grid_points[:size]

        return grid_points

    def write_param_file(self):
        """
        Writes the sampled parameters to a file specified in the configuration.
        Supports CSV and FITS file formats.
        """
        directory = self.config["sampling"]["output_directory"]
        filename = directory + self.config["sampling"]["output_file"]
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
        """
        Fills in missing m_v, surface_brightness, or r_scale based on other parameters.

        Parameters:
        - calc_param: The parameter to calculate.
        - df: DataFrame with the dwarf galaxies' parameters.

        Returns:
        - An array of calculated parameter values.
        """
        if calc_param == "m_v":
            vals = sb_rh_to_mv(
                sb=df["surface_brightness"],
                rh=df["r_scale"],
                dist=df["dist"],
            )
        elif calc_param == "r_scale":
            vals = sb_mv_to_rh(
                sb=df["surface_brightness"], M_v=df["m_v"], dist=df["dist"]
            )
        elif calc_param == "surface_brightness":
            vals = rh_mv_to_sb(M_v=df["m_v"], rh=df["r_scale"], dist=df["dist"])
        else:
            raise Exception(f"unknown calc param: {calc_param}")
        return vals

    def _populate_ra_dec(self):
        """
        Populates the RA and Dec fields of the dwarf galaxies based on their x/y_cen.

        Returns:
        - Arrays of RA and Dec values.
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
