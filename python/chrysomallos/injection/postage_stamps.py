import os
import time

import fitsio
import lsst.source.injection as si
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import make_lupton_rgb
from tqdm import tqdm

from chrysomallos.utils import logger

__all__ = [
    "PostageStampGenerator",
]


class PostageStampGenerator:
    """
    injectes dwarf galaxies into images and generates postage stamps of.
    """

    def __init__(self, config, dwarf_params_frame, dwarf_catalogs, coadd_dict) -> None:
        """
        Initializes the generator with configurations and necessary data structures.

        Parameters:
        - config: Configuration dictionary for stamp generation parameters.
        - dwarf_params_frame: DataFrame containing parameters of dwarf galaxies.
        - dwarf_catalogs: Dictionary of catalogs of dwarf galaxies.
        - coadd_dict: Dictionary containing coadded image data.
        """
        self.config = config
        self.dwarf_catalogs = dwarf_catalogs
        self.coadd_dict = coadd_dict
        self.dwarf_params_frame = dwarf_params_frame

        # check that directory exists
        os.makedirs(self.config["stamp"]["directory"], exist_ok=True)

        # make sure bands is a list with 3 elements
        bands = self.config["pipelines"]["bands"]
        if len(self.config["pipelines"]["bands"]) != 3:
            raise ValueError(
                f"3 bands must be specified for postage stamp generation : {bands}"
            )

    def run(self):
        """
        Executes the process of generating and saving postage stamps
        for each dwarf galaxy in the input frame.
        """
        
        # setup injection task
        inject_config = si.CoaddInjectConfig()
        inject_task = si.CoaddInjectTask(config=inject_config)

        stamp_x_size = self.config["stamp"]["size"][0]
        stamp_y_size = self.config["stamp"]["size"][1]

        first_band = self.config["pipelines"]["bands"][0]

        fig, ax = self.create_axes_for_stamps(
            stamp_x_size=stamp_x_size, stamp_y_size=stamp_y_size, dpi=100
        )

        for i in tqdm(range(self.config["sampling"]["n_dwarfs"])):
            injection_dict = {}
            # number of sources we are injecting
            cat_length = len(self.dwarf_catalogs[first_band][i])

            x_offset = np.random.uniform(-int(stamp_x_size / 2.2), int(stamp_x_size / 2.2))
            y_offset = np.random.uniform(-int(stamp_y_size / 2.2), int(stamp_y_size / 2.2))
            
            # make title of stamp
            title = self.make_stamp_title(
                stamp_directory=self.config["stamp"]["directory"],
                stamp_title_prefix=self.config["stamp"]["title_format"],
                tract=self.config["pipelines"]["tract"],
                patch=self.config["pipelines"]["patch"],
                dwarf_id=i,
                ra=self.dwarf_params_frame["ra"][i],
                dec=self.dwarf_params_frame["dec"][i],
                bands=self.config["pipelines"]["bands"],
                x_offset=x_offset,
                y_offset=y_offset,
            )

            logger.info(f"creating {title} with {cat_length} sources.")
            
            
            
            start_time = time.time()
            for band in self.config["pipelines"]["bands"]:
                image = self.coadd_dict[band]["image"]
                input_exposure = image.clone()
                psf = self.coadd_dict[band]["psf"]
                photo_calib = self.coadd_dict[band]["photo_calib"]
                wcs = self.coadd_dict[band]["wcs"]
                bbox = self.coadd_dict[band]["bbox"]

                x_cen = self.dwarf_params_frame["x_cen"][i]
                y_cen = self.dwarf_params_frame["y_cen"][i]

                minx = bbox.beginX
                miny = bbox.beginY

                # crop injection catalog to stamp size
                injection_catalogs = self.crop_injection_catalog(
                    catalog=self.dwarf_catalogs[band][i],
                    band=band,
                    x_cen=x_cen,
                    y_cen=y_cen,
                    stamp_x_size=stamp_x_size,
                    stamp_y_size=stamp_y_size,
                )

                injection_catalogs = [injection_catalogs]
                inject_output = inject_task.run(
                    injection_catalogs=injection_catalogs,
                    input_exposure=input_exposure,
                    psf=psf,
                    photo_calib=photo_calib,
                    wcs=wcs,
                )
                exposure = inject_output.output_exposure
                
                stamp_range = [
                    minx + x_cen - int(stamp_x_size / 2) + x_offset,
                    minx + x_cen + int(stamp_x_size / 2) + x_offset,
                    miny + y_cen - int(stamp_y_size / 2) + y_offset,
                    miny + y_cen + int(stamp_y_size / 2) + y_offset
                ]

                injection_dict[band] = exposure.image[
                    stamp_range[0] : stamp_range[1], stamp_range[2] : stamp_range[3]
                ].array

            self.make_one_stamp_png(
                injection_dict=injection_dict,
                title=title,
                Q=self.config["stamp"]["Q"],
                stretch=self.config["stamp"]["stretch"],
                minimum=self.config["stamp"]["minimum"],
                bands=self.config["pipelines"]["bands"],
                ax=ax,
            )
            end_time = time.time()
            logger.info(
                f"Time to create stamp with {cat_length} sources: {end_time - start_time:0.2f} seconds"
            )

    def crop_injection_catalog(
        self, catalog, band, x_cen, y_cen, stamp_x_size, stamp_y_size
    ):
        """
        Crops the injection catalog to the size of the postage stamp.

        Parameters:
        - catalog: The catalog of sources to be injected.
        - band: The band of observation.
        - x_cen, y_cen: Central coordinates of the postage stamp.
        - stamp_x_size, stamp_y_size: Dimensions of the postage stamp.

        Returns:
        - Cropped catalog of sources.
        """
        wcs = self.coadd_dict[band]["wcs"]
        bbox = self.coadd_dict[band]["bbox"]

        x_pix, y_pix = wcs.skyToPixelArray(catalog["ra"], catalog["dec"], degrees=True)
        sel = (
            (abs(x_pix - x_cen - bbox.beginX) < stamp_x_size / 2)
            & (abs(y_pix - y_cen - bbox.beginY) < stamp_y_size)
        ) | (catalog["source_type"] == "Sersic")
        cropped_catalog = catalog[sel]
        return cropped_catalog

    def make_one_stamp_png(
        self, injection_dict, title, Q, stretch, minimum, bands, ax=None
    ):
        """
        Generates a single postage stamp image and saves it as a PNG file.

        Parameters:
        - injection_dict: Dictionary of image arrays by band.
        - title: Filename for the saved PNG image.
        - Q, stretch, minimum: Parameters for Lupton RGB image creation.
        - ax: Matplotlib Axes object for plotting (optional).
        """
        if ax is None:
            fig, ax = self.create_axes_for_stamps()
        else:
            ax.clear()
        rgb = make_lupton_rgb(
            injection_dict[bands[2]],
            injection_dict[bands[1]],
            injection_dict[bands[0]],
            Q=Q,
            stretch=stretch,
            minimum=minimum,
        )
        ax.imshow(rgb)

        plt.savefig(title, dpi=100)

    def create_axes_for_stamps(self, stamp_x_size=600, stamp_y_size=600, dpi=100):
        """
        Creates Matplotlib axes suitable for postage stamp images.

        Parameters:
        - stamp_x_size, stamp_y_size: Dimensions of the postage stamp.
        - dpi: Dots per inch for the output image.

        Returns:
        - fig, ax: Matplotlib figure and axes objects.
        """
        x_size_inches = int(stamp_x_size / dpi)
        y_size_inches = int(stamp_y_size / dpi)
        fig, ax = plt.subplots(figsize=(x_size_inches, y_size_inches), dpi=dpi)
        ax.axis("off")
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        ax.set_xticks([], [])
        ax.set_yticks([], [])
        return fig, ax

    def make_stamp_title(
        self,
        stamp_directory,
        stamp_title_prefix,
        tract,
        patch,
        dwarf_id,
        ra,
        dec,
        bands,
        x_offset,
        y_offset,
    ):
        """
        Generates a title for a postage stamp image file.

        Parameters:
        - stamp_directory: Directory where the stamp will be saved.
        - stamp_title_prefix: Prefix for the stamp's filename.
        - tract, patch: Tract and patch identifiers.
        - dwarf_id: Identifier for the dwarf galaxy.
        - ra, dec: Right ascension and declination of the galaxy.

        Returns:
        - Filename for the stamp image.
        """
        band_str = "".join(bands)
        filename = stamp_directory + stamp_title_prefix + "_"
        filename += f"{tract}_{patch}_{dwarf_id}_{ra:0.2f}_{dec:0.2f}_{band_str}_offset_{x_offset:0.0f}_{y_offset:0.0f}.png"
        return filename

    def save_stamp_as_fits(self, stamp_directory, title, injection_dict):
        """
        Saves the generated postage stamp as a FITS file.

        Parameters:
        - stamp_directory: Directory to save the FITS file.
        - title: Title of the stamp, used for the FITS filename.
        - injection_dict: Dictionary of image data by band.
        """
        fits = fitsio.FITS(stamp_directory + title.replace(".png", ".fits"), "rw")
        array_list = [
            injection_dict[band] for band in self.config["pipelines"]["bands"]
        ]
        names = [band for band in self.config["pipelines"]["bands"]]
        fits.write(array_list, names=names, clobber=True)
        fits.close()

    def generate_empty_stamps(self, n_stamps):
        """
        Generates a specified number of empty postage stamps.

        Parameters:
        - n_stamps: Number of empty stamps to generate.
        """
        stamp_x_size = self.config["stamp"]["size"][0]
        stamp_y_size = self.config["stamp"]["size"][1]
        tract = self.config["pipelines"]["tract"]
        patch = self.config["pipelines"]["patch"]
        first_band = self.config["pipelines"]["bands"][0]
        wcs = self.coadd_dict[first_band]["wcs"]

        fig, ax = self.create_axes_for_stamps(
            stamp_x_size=stamp_x_size, stamp_y_size=stamp_y_size, dpi=100
        )
        x_cen = np.random.uniform(
            self.config["sampling"]["params"]["x_cen"][0],
            self.config["sampling"]["params"]["x_cen"][1],
            n_stamps,
        )
        y_cen = np.random.uniform(
            self.config["sampling"]["params"]["y_cen"][0],
            self.config["sampling"]["params"]["y_cen"][1],
            n_stamps,
        )
        # get ra/dec from wcs
        ra, dec = wcs.pixelToSkyArray(x_cen, y_cen, degrees=True)
        for i in range(n_stamps):
            injection_dict = {}

            for band in self.config["pipelines"]["bands"]:
                # make title
                title = (
                    self.config["stamp"]["directory"]
                    + f"empty_stamp_{tract}_{patch}_{i}_{ra[i]:0.2f}_{dec[i]:0.2f}.png"
                )

                # stamp parameters
                minx = self.coadd_dict[band]["bbox"].beginX
                miny = self.coadd_dict[band]["bbox"].beginY
                stamp_range = [
                    minx + x_cen[i] - int(stamp_x_size / 2),
                    minx + x_cen[i] + int(stamp_x_size / 2),
                    miny + y_cen[i] - int(stamp_y_size / 2),
                    miny + y_cen[i] + int(stamp_y_size / 2),
                ]

                injection_dict[band] = (
                    self.coadd_dict[band]["image"]
                    .clone()
                    .image[
                        stamp_range[0] : stamp_range[1], stamp_range[2] : stamp_range[3]
                    ]
                    .array
                )

            self.make_one_stamp_png(
                injection_dict=injection_dict,
                title=title,
                Q=self.config["stamp"]["Q"],
                stretch=self.config["stamp"]["stretch"],
                minimum=self.config["stamp"]["minimum"],
                bands=self.config["pipelines"]["bands"],
                ax=ax,
            )
