import astropy.table as atable
import lsst.source.injection as si
import numpy as np
import yaml
from lsst.daf.butler import Butler, CollectionType
from lsst.daf.butler.registry import MissingCollectionError

from chrysomallos.injection import adopt_a_cat, massage_the_cat
from chrysomallos.utils import get_coadd_dict, logger


class CreateDwarfInjectionCatalog:
    """
    Class to create an injection catalog for dwarf galaxies.
    Processes configuration data, fetches coadd from butler,
    and creates catalogs ready for injection, then ingests into butler.
    """

    def __init__(self, config, dwarf_params_frame, coadd_dict=None):
        """
        Initialize the class with the given configuration dictionary.

        Parameters
        ----------
        config_dict : dict, optional
            Configuration dictionary with parameters and settings.
            Default is `default_config_dict`.
        """
        self.config = config
        self.dwarf_params_frame = dwarf_params_frame.reset_index(
            drop=True, inplace=False
        )

        self.coadd_dict = get_coadd_dict(coadd_dict, config)
        self.dwarf_cats = []
        self.first_band = config["pipelines"]["bands"][0]

        if (len(np.unique(dwarf_params_frame["tract"]))) > 1:
            raise Exception(
                f"can only run one tract at at time: {np.unique(dwarf_params_frame['tract'])}"
            )
        if (len(np.unique(dwarf_params_frame["patch"]))) > 1:
            raise Exception(
                f"can only run one patch at at time: {np.unique(dwarf_params_frame['patch'])}"
            )

    def run(self, ingest=False, multiproc=False):
        """
        Main method to run the dwarf injection catalog creation process.
        """
        # generate xy catalog for each dwarf

        self.generate_catalogs(multiproc=multiproc)

        logger.info(f"generated catalogs for {len(self.dwarf_params_frame)} dwarfs")

        if self.config["injection"]["type"] == "grid":
            # convert to ssi format and concatenate to single catalog per band
            self.injection_cats = {}
            for band in self.config["pipelines"]["bands"]:
                self.injection_cats[band] = []
            for i, cat in enumerate(self.dwarf_cats):
                for band in self.config["pipelines"]["bands"]:
                    self.injection_cats[band].append(
                        massage_the_cat(
                            cat_inp=cat,
                            mag_limit=self.config["injection"]["mag_limit"],
                            band_for_injection=band,
                            wcs=self.coadd_dict[self.first_band]["wcs"],
                            bbox=self.coadd_dict[self.first_band]["bbox"],
                            x_cen=self.dwarf_params_frame["x_cen"][i],
                            y_cen=self.dwarf_params_frame["y_cen"][i],
                            r_scale=self.dwarf_params_frame["r_scale"][i],
                            dist=self.dwarf_params_frame["distance"][i],
                        )
                    )
            for band in self.config["pipelines"]["bands"]:
                self.injection_cats[band] = atable.vstack(self.injection_cats[band])
            if ingest:
                self.ingest_injection_catalogs(
                    si_input_collection=self.inject_cat_collection,
                    catalogs=self.injection_cats,
                    bands=self.config["pipelines"]["bands"],
                )
            else:
                logger.info(
                    "not ingesting catalogs to change use self.run(ingest=True)"
                )

        elif self.config["injection"]["type"] == "stamp":
            # for each catalog in self.dwarf_cats, run massage_the_cat and
            # save result as a dict where is key is the config never ingest
            self.injection_cats = {
                band: {} for band in self.config["pipelines"]["bands"]
            }

            for i, cat in enumerate(self.dwarf_cats):
                for band in self.config["pipelines"]["bands"]:
                    self.injection_cats[band][i] = massage_the_cat(
                        cat_inp=cat,
                        mag_limit=self.config["injection"]["mag_limit"],
                        band_for_injection=band,
                        wcs=self.coadd_dict[band]["wcs"],
                        bbox=self.coadd_dict[band]["bbox"],
                        x_cen=self.dwarf_params_frame["x_cen"][i],
                        y_cen=self.dwarf_params_frame["y_cen"][i],
                        r_scale=self.dwarf_params_frame["r_scale"][i],
                        dist=self.dwarf_params_frame["distance"][i],
                    )
        else:
            raise Exception(
                f"unknown injection type: {self.config['injection']['type']}"
            )

        return self.injection_cats, self.coadd_dict

    def generate_catalogs(self, multiproc=False):
        self.config["injection"]["mag_limit_band"]

        rows = self.dwarf_params_frame.to_dict(orient="records")
        r_len = len(rows)
        gen_args = list(
            zip(
                rows,
                r_len * [self.coadd_dict[self.config["pipelines"]["bands"][0]]["wcs"]],
                r_len * [self.coadd_dict[self.config["pipelines"]["bands"][0]]["bbox"]],
                r_len * [self.config["injection"]["mag_limit_band"]],
            )
        )

        if multiproc:
            from multiprocessing import Pool

            processes = multiproc if multiproc > 0 else None
            p = Pool(processes, maxtasksperchild=1)
            out = p.map(self.generate_single_catalog, gen_args)
        else:
            out = [self.generate_single_catalog(args) for args in gen_args]

        self.dwarf_cats = out

    @staticmethod
    def generate_single_catalog(args):
        row, wcs, bbox, mag_limit_band = args

        if np.isnan(row["random_seed_injection"]) | (
            row["random_seed_injection"] is None
        ):
            random_seed = None
        else:
            random_seed = int(row["random_seed_injection"])

        catalog = adopt_a_cat(
            wcs=wcs,
            bbox=bbox,
            age=row["age"],
            feh=row["feh"],
            stellar_mass=row["stellar_mass"],
            dist=row["distance"],
            r_scale=row["r_scale"],
            ellip=row["ellipticity"],
            theta=row["theta"],
            n=row["n"],
            m_v=row["m_v"],
            mag_limit=36,  # notice this is set to be deeper
            mag_limit_band=mag_limit_band,
            random_seed=random_seed,
        )
        return catalog

    # def get_data_ids(self, tract=9615, patch=3):
    #     """
    #     Get data IDs based on band, skymap, tract, and patch information.

    #     Parameters
    #     ----------
    #     tract : int, optional
    #         Tract number. Default is 9615.
    #     patch : int, optional
    #         Patch number. Default is 3.
    #     """
    #     self.dataid_dict = {}
    #     for band in self.config["pipelines"]["bands"]:
    #         self.dataid_dict[band] = {
    #             "band": band,
    #             "skymap": "hsc_rings_v1",
    #             "tract": tract,
    #             "patch": patch,
    #         }

    # def get_coadds(self):
    #     """
    #     Fetch coadd data based on the specified bands and data IDs.
    #     """
    #     # TODO: think about converting to static method and putting in utils
    #     self.coadd_dict = {band: {} for band in self.config["pipelines"]["bands"]}
    #     for band in self.config["pipelines"]["bands"]:
    #         image = self.butler.get("deepCoadd_calexp", dataId=self.dataid_dict[band])
    #         self.coadd_dict[band]["image"] = image
    #         self.coadd_dict[band]["wcs"] = image.getWcs()
    #         self.coadd_dict[band]["bbox"] = image.getBBox()
    #         self.coadd_dict[band]["psf"] = image.getPsf()
    #         self.coadd_dict[band]["photo_calib"] = image.getPhotoCalib()

    def ingest_injection_catalogs(self, si_input_collection, catalogs, bands):
        """
        Ingest the prepared injection catalogs into the specified collection.

        Parameters
        ----------
        si_input_collection : str
            The collection to ingest the catalogs into.
        catalogs : dict
            Dictionary of prepared catalogs, keyed by band.
        bands : list
            List of bands to consider.
        """
        writeable_butler = Butler(self.repo, writeable=True)
        try:
            writeable_butler.removeRuns([si_input_collection])
        except MissingCollectionError:
            logger.info("Writing into a new RUN collection")
            pass
        else:
            logger.info("Prior RUN collection located and successfully removed")
        _ = writeable_butler.registry.registerCollection(
            si_input_collection, type=CollectionType.RUN
        )

        self.injection_cat_refs = {}
        for band in bands:
            self.injection_cat_refs[band] = si.ingest_injection_catalog(
                writeable_butler=writeable_butler,
                table=catalogs[band],
                band=band,
                output_collection=si_input_collection,
                dataset_type_name="injection_catalog",
            )

    def save_config(self, filename):
        yaml_dict = vars(self).copy()
        del yaml_dict["dwarf_configs"]
        with open(filename, "w") as file:
            yaml.dump(yaml_dict, file)


# if __name__ == "__main__":
#     """
#     Main execution script.
#     Parses input arguments, reads the configuration file,
#     and runs the dwarf injection catalog creation.
#     """
#     config_dict = default_config_dict.copy()
#     parser = argparse.ArgumentParser(
#         description="Description of your script/functionality."
#     )
#     parser.add_argument("filename", help="name of config.yaml file")
#     args = parser.parse_args()

#     with open(args.filename, "r") as file:
#         user_config_dict = yaml.safe_load(file)
#         for key, value in user_config_dict.items():
#             config_dict[key] = value

#     creator = CreateDwarfInjectionCatalog(config_dict)
#     creator.run(ingest=False)
