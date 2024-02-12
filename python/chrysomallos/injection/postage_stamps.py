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
    def __init__(self, config, dwarf_params_frame, dwarf_catalogs, coadd_dict) -> None:
        self.config = config
        self.dwarf_catalogs = dwarf_catalogs
        self.coadd_dict = coadd_dict
        self.dwarf_params_frame = dwarf_params_frame

    def run(self):
        # check that directory exists
        os.makedirs(self.config["stamp"]["directory"], exist_ok=True)

        # setup injection task
        inject_config = si.CoaddInjectConfig()
        inject_task = si.CoaddInjectTask(config=inject_config)

        stamp_x_size = self.config["stamp"]["size"][0]
        stamp_y_size = self.config["stamp"]["size"][1]

        fig, ax = self.create_axes_for_stamps(
            stamp_x_size=stamp_x_size, stamp_y_size=stamp_y_size, dpi=100
        )

        for i in tqdm(range(self.config["sampling"]["n_dwarfs"])):
            injection_dict = {}
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
                    minx + x_cen - int(stamp_x_size / 2),
                    minx + x_cen + int(stamp_x_size / 2),
                    miny + y_cen - int(stamp_y_size / 2),
                    miny + y_cen + int(stamp_y_size / 2),
                ]

                injection_dict[band] = exposure.image[
                    stamp_range[0] : stamp_range[1], stamp_range[2] : stamp_range[3]
                ].array

            title = self.make_stamp_title(
                stamp_directory=self.config["stamp"]["directory"],
                stamp_title_prefix=self.config["stamp"]["title_format"],
                tract=self.config["pipelines"]["tract"],
                patch=self.config["pipelines"]["patch"],
                dwarf_id=i,
                ra=self.dwarf_params_frame["ra"][i],
                dec=self.dwarf_params_frame["dec"][i],
            )

            self.make_one_stamp_png(
                injection_dict=injection_dict,
                title=title,
                Q=self.config["stamp"]["Q"],
                stretch=self.config["stamp"]["stretch"],
                minimum=self.config["stamp"]["minimum"],
                ax=ax,
            )
            end_time = time.time()
            cat_length = len(self.dwarf_catalogs[band][i])
            logger.info(
                f"Time to create stamp with {cat_length} sources: {end_time - start_time:0.2f} seconds"
            )

    def generate_empty_stamps(self, n_stamps):
        stamp_x_size = self.config["stamp"]["size"][0]
        stamp_y_size = self.config["stamp"]["size"][1]
        tract = self.config["pipelines"]["tract"]
        patch = self.config["pipelines"]["patch"]

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
        for i in range(n_stamps):
            injection_dict = {}

            for band in self.config["pipelines"]["bands"]:
                # get ra/dec from wcs
                ra, dec = self.coadd_dict[band]["wcs"].pixelToSky(
                    float(x_cen[i]), float(y_cen[i])
                )
                # make title
                title = (
                    self.config["stamp"]["directory"]
                    + f"empty_stamp_{tract}_{patch}_{i}_{ra}_{dec}.png"
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
                ax=ax,
            )

    def crop_injection_catalog(
        self, catalog, band, x_cen, y_cen, stamp_x_size, stamp_y_size
    ):
        wcs = self.coadd_dict[band]["wcs"]
        bbox = self.coadd_dict[band]["bbox"]

        x_pix, y_pix = wcs.skyToPixelArray(catalog["ra"], catalog["dec"], degrees=True)
        sel = (
            (abs(x_pix - x_cen - bbox.beginX) < stamp_x_size / 2)
            & (abs(y_pix - y_cen - bbox.beginY) < stamp_y_size)
        ) | (catalog["source_type"] == "Sersic")
        cropped_catalog = catalog[sel]
        return cropped_catalog

    def make_one_stamp_png(self, injection_dict, title, Q, stretch, minimum, ax=None):
        if ax is None:
            fig, ax = self.create_axes_for_stamps()
        else:
            ax.clear()
        rgb = make_lupton_rgb(
            injection_dict["i"],
            injection_dict["r"],
            injection_dict["g"],
            Q=Q,
            stretch=stretch,
            minimum=minimum,
        )
        ax.imshow(rgb)

        plt.savefig(title, dpi=100)

    def create_axes_for_stamps(self, stamp_x_size=600, stamp_y_size=600, dpi=100):
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
    ):
        filename = stamp_directory + stamp_title_prefix
        filename += f"{tract}_{patch}_{dwarf_id}_{ra:0.2f}_{dec:0.2f}.png"
        return filename

    def save_stamp_as_fits(self, stamp_directory, title, injection_dict):
        fits = fitsio.FITS(stamp_directory + title.replace(".png", ".fits"), "rw")
        array_list = [
            injection_dict[band] for band in self.config["pipelines"]["bands"]
        ]
        names = [band for band in self.config["pipelines"]["bands"]]
        fits.write(array_list, names=names, clobber=True)
        fits.close()


# class CatalogCreator:
#     def __init__(
#         self, butler, bands=["g"], default_config_dict=None, default_dwarf_dict=None
#     ):
#         self.butler = butler
#         self.bands = bands
#         self.default_config_dict = default_config_dict or {}
#         self.default_dwarf_dict = default_dwarf_dict or {}
#         self.catalogs = {}
#         self.image_dict = self._initialize_image_dict()

#     def _initialize_image_dict(self):
#         image_dict = {}
#         for band in self.bands:
#             dataid = {"band": band, "skymap": "hsc_rings_v1", "tract": 9615, "patch": 3}
#             image_dict[band] = self.butler.get("deepCoadd_calexp", dataId=dataid)
#         return image_dict

#     def create_catalogs(
#         self, surface_brighness_vals, mv_vals, dist, x_cen, y_cen, mag_lim
#     ):
#         wcs = self.image_dict["g"].getWcs()
#         bbox = self.image_dict["g"].getBBox()
#         j = 0

#         for sb in tqdm(surface_brighness_vals[:]):
#             for m_v in mv_vals:
#                 r_h = sb_mv_to_rh(
#                     sb, m_v, distance=2e6
#                 )  # Assuming this function is defined elsewhere
#                 stellar_mass = mstar_from_absmag(
#                     m_v
#                 )  # Assuming this function is defined elsewhere

#                 dwarf_dicts = []
#                 new_config_dict = self.default_config_dict.copy()
#                 new_config_dict["inject_cat_collection"] = "test"

#                 dwarf_dict = self.default_dwarf_dict.copy()
#                 dwarf_dict.update(
#                     {
#                         "id": j,
#                         "dist": dist,
#                         "x_cen": x_cen,
#                         "y_cen": y_cen,
#                         "r_scale": r_h,
#                         "sb": sb,
#                         "mag_limit": mag_lim,
#                         "m_v": m_v,
#                         "stellar_mass": stellar_mass,
#                     }
#                 )

#                 dwarf_dicts.append(dwarf_dict)
#                 new_config_dict["dwarfs"] = dwarf_dicts

#                 creator = CreateDwarfInjectionCatalog(new_config_dict)
#                 tmp_catalog = creator.run(ingest=False, coadd_dict=self.image_dict)["g"]

#                 x_pix, y_pix = wcs.skyToPixelArray(
#                     tmp_catalog["ra"], tmp_catalog["dec"], degrees=True
#                 )
#                 sel = (
#                     (abs(x_pix - x_cen - bbox.beginX) < 400)
#                     & (abs(y_pix - y_cen - bbox.beginY) < 400)
#                 ) | (tmp_catalog["source_type"] == "Sersic")
#                 new_catalog_dict = {
#                     "config": new_config_dict,
#                     "catalog": tmp_catalog[sel],
#                 }

#                 self.catalogs[j] = new_catalog_dict
#                 j += 1

#     def get_catalogs(self):
#         return self.catalogs

#     def _create_catalog_for_params(self, params):
#         sb, m_v, wcs, bbox, dist, x_cen, y_cen, mag_lim = params

#         r_h = sb_mv_to_rh(sb, m_v, distance=2e6)
#         stellar_mass = mstar_from_absmag(m_v)

#         dwarf_dicts = []
#         new_config_dict = self.default_config_dict.copy()
#         new_config_dict["inject_cat_collection"] = "test"

#         dwarf_dict = self.default_dwarf_dict.copy()
#         dwarf_dict.update(
#             {
#                 "id": self.j,
#                 "dist": dist,
#                 "x_cen": x_cen,
#                 "y_cen": y_cen,
#                 "r_scale": r_h,
#                 "sb": sb,
#                 "mag_limit": mag_lim,
#                 "m_v": m_v,
#                 "stellar_mass": stellar_mass,
#             }
#         )

#         dwarf_dicts.append(dwarf_dict)
#         new_config_dict["dwarfs"] = dwarf_dicts

#         creator = CreateDwarfInjectionCatalog(new_config_dict)
#         tmp_catalog = creator.run(ingest=False, coadd_dict=self.image_dict)["g"]

#         x_pix, y_pix = wcs.skyToPixelArray(
#             tmp_catalog["ra"], tmp_catalog["dec"], degrees=True
#         )
#         sel = (
#             (abs(x_pix - x_cen - bbox.beginX) < 400)
#             & (abs(y_pix - y_cen - bbox.beginY) < 400)
#         ) | (tmp_catalog["source_type"] == "Sersic")
#         new_catalog_dict = {"config": new_config_dict, "catalog": tmp_catalog[sel]}

#         return (self.j, new_catalog_dict)

#     # def create_catalogs(
#     #     self, surface_brighness_vals, mV_vals, dist, x_cen, y_cen, mag_lim
#     # ):
#     #     wcs = self.image_dict["g"].getWcs()
#     #     bbox = self.image_dict["g"].getBBox()

#     #     params = [
#     #         (sb, mV, wcs, bbox, dist, x_cen, y_cen, mag_lim)
#     #         for sb in surface_brighness_vals
#     #         for mV in mV_vals
#     #     ]

#     #     with multiprocessing.Pool() as pool:
#     #         results = list(
#     #             tqdm(
#     #                 pool.imap(self._create_catalog_for_params, params),
#     #                 total=len(params),
#     #             )
#     #         )

#     #     for j, catalog_dict in results:
#     #         self.catalogs[j] = catalog_dict
#     #         self.j += 1
