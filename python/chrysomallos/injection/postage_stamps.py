# import logging

import fitsio
import lsst.source.injection as si
from tqdm import tqdm

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
        # for each catalog in dwarf catalogs inject in each band and save fits file with each band as
        # a different extension use fitsio to save the fits file
        inject_config = si.CoaddInjectConfig()
        inject_task = si.CoaddInjectTask(config=inject_config)

        for i in tqdm(range(self.config["sampling"]["n_dwarfs"])):
            injection_dict = {}
            for band in self.config["pipelines"]["bands"]:
                image = self.coadd_dict[band]["image"]
                input_exposure = image.clone()
                psf = image.getPsf()
                photo_calib = image.getPhotoCalib()
                wcs = self.coadd_dict[band]["wcs"]
                bbox = self.coadd_dict[band]["bbox"]

                x_cen = self.dwarf_params_frame["x_cen"][i]
                y_cen = self.dwarf_params_frame["y_cen"][i]
                stamp_x_size = self.config["injection"]["stamp_size"][0]
                stamp_y_size = self.config["injection"]["stamp_size"][1]

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

            fits = fitsio.FITS("test.fits", "rw")
            array_list = [
                injection_dict[band] for band in self.config["pipelines"]["bands"]
            ]
            names = [band for band in self.config["pipelines"]["bands"]]
            fits.write(array_list, names=names, clobber=True)

            # self.dwarf_catalogs[i]["image"] = inject_output.output_exposure

        # injection
        # for i in tqdm(self.catalogs):
        #     inject_output = inject_task.run(
        #         injection_catalogs=[catalogs[i]["catalog"]],
        #         input_exposure=coadd_dict["g"].clone(),
        #         psf=coadd_dict["g"].getPsf(),
        #         photo_calib=coadd_dict["g"].getPhotoCalib(),
        #         wcs=coadd_dict["g"].getWcs(),
        #     )
        #     catalogs[i]["image"] = inject_output.output_exposure

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
