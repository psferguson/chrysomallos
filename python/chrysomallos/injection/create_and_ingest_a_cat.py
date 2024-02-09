import argparse

import astropy.table as atable
import lsst.source.injection as si
import numpy as np
import yaml
from lsst.daf.butler import Butler, CollectionType
from lsst.daf.butler.registry import MissingCollectionError

from chrysomallos.injection import adopt_a_cat, massage_the_cat
from chrysomallos.utils.log import logger

# from dataclasses import dataclass
# from chrysomallos.utils.utils import (
#     mstar_from_absmag,
#     rh_mv_to_sb,
#     sb_mv_to_rh,
#     sb_rh_to_mv,
# )

default_config_dict = {
    "repo": "/repo/main",
    "collection": "HSC/runs/RC2/w_2023_32/DM-40356",
    "bands": ["g", "r", "i"],
    "tract": 9615,
    "patch": 3,
    "dwarfs": [],
}


# @dataclass
# class DwarfConfig:
#     """
#     Configuration class for a dwarf galaxy.
#     Sets various properties of the dwarf galaxy, including its location,
#     age, metallicity, and other parameters.
#     """

#     # DEFAULT PARAMETERS
#     # ID of the injection cat
#     id: int
#     # X-coordinate of the object's center [pixel]
#     x_cen: float
#     # Y-coordinate of the object's center [pixel]
#     y_cen: float
#     # surface brightness
#     sb: float = np.nan
#     # absolute magnitude
#     m_v: float = np.nan
#     # Age of the object [Gyr]
#     age: float = 10.0
#     # Metallicity of the object [Fe/H]
#     feh: float = -2.0
#     # Stellar mass of the object [Msolar]
#     stellar_mass: float = 5.0e5
#     # Distance to the object [Mpc]
#     dist: float = np.nan
#     # Scale radius of the object [pc]
#     r_scale: float = np.nan
#     # Ellipticity of the object (0-1)
#     ellip: float = 0
#     # Orientation angle of the object
#     theta: float = 0
#     # Sersic index of the object (shape parameter)
#     n: float = 1
#     # Magnitude limit to inject for the object
#     mag_limit: float = 27
#     # Band in which the magnitude limit is defined
#     mag_limit_band: str = "LSST_g"
#     # Seed for random number generation (for reproducibility)
#     random_seed: int = None
#     #

#     def __post_init__(self, dwarf_dict=None, **kwargs):
#         """
#         Initialize the class with the given dwarf dictionary.

#         Parameters
#         ----------
#         dwarf_dict : dict
#             Dictionary containing the properties and attributes of a dwarf.
#         """
#         # check that values are not bad
#         if (self.ellip < 0) or (self.ellip > 1):
#             raise RuntimeError(f"ellip must be between 0 and 1 : {self.ellip}")
#         if (
#             np.isnan(self.Mv)
#             and not np.isnan(self.r_scale)
#             and not np.isnan(self.sb)
#             and not np.isnan(self.dist)
#         ):
#             self.Mv = sb_rh_to_mv(sb=self.sb, rh=self.r_scale, distance=self.dist)
#         if (
#             np.isnan(self.r_scale)
#             and not np.isnan(self.Mv)
#             and not np.isnan(self.sb)
#             and not np.isnan(self.dist)
#         ):
#             self.Mv = sb_mv_to_rh(sb=self.sb, M_v=self.Mv, distance=self.dist)
#         if (
#             np.isnan(self.sb)
#             and not np.isnan(self.Mv)
#             and not np.isnan(self.r_scale)
#             and not np.isnan(self.dist)
#         ):
#             self.sb = rh_mv_to_sb(rh=self.r_scale, M_v=self.Mv, distance=self.dist)
#         if np.isnan(self.mass) and not np.isnan(self.Mv):
#             self.mass = mstar_from_absmag(self.mV)
#     @classmethod
#     def from_config_and_frame(self, config, dwarf_params_frame):
#         # fetch the constructor's signature
#         cls_fields = {field for field in signature(cls).parameters}

#         for field in cls_fields:

#         # split the kwargs into native ones and new ones
#         native_args, new_args = {}, {}
#         for name, val in kwargs.items():
#             if name in cls_fields:
#                 native_args[name] = val
#             else:
#                 new_args[name] = val

#         # use the native ones to create the class ...
#         ret = cls(**native_args)

#         # ... and add the new ones by hand
#         for new_name, new_val in new_args.items():
#             setattr(ret, new_name, new_val)
#         return ret


class CreateDwarfInjectionCatalog:
    """
    Class to create an injection catalog for dwarf galaxies.
    Processes configuration data, fetches coadd from butler,
    and creates catalogs ready for injection, then ingests into butler.
    """

    def __init__(self, config, dwarf_params_frame):
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

        if (len(np.unique(dwarf_params_frame["tract"]))) > 1:
            raise Exception(
                f"can only run one tract at at time: {np.unique(dwarf_params_frame['tract'])}"
            )
        if (len(np.unique(dwarf_params_frame["patch"]))) > 1:
            raise Exception(
                f"can only run one patch at at time: {np.unique(dwarf_params_frame['patch'])}"
            )
        # for key, value in config_dict.items():
        #     setattr(self, key, value)
        # self.dwarf_configs = []  # this should be a dataframe
        # for dwarf in self.dwarfs:
        #     self.dwarf_configs.append(DwarfConfig(**dwarf))

    def run(self, ingest=False, coadd_dict=None, multiproc=False):
        """
        Main method to run the dwarf injection catalog creation process.
        """

        self.get_data_ids(
            tract=self.config["pipelines"]["tract"],
            patch=self.config["pipelines"]["patch"],
        )

        if coadd_dict:
            first_band = self.config["pipelines"]["bands"][0]
            self.coadd_dict = coadd_dict
            self.wcs = self.coadd_dict[first_band].getWcs()
            self.bbox = self.coadd_dict[first_band].getBBox()
        else:
            self.butler = Butler(
                self.config["pipelines"]["repo"],
                collections=self.config["pipelines"]["input_collections"],
            )
            # grab deepCoaddd
            self.get_coadds()
        self.dwarf_cats = []
        # generate xy catalog for each dwarf

        self.generate_catalogs(multiproc=multiproc)

        logger.info(f"generated catalogs for {len(self.dwarf_params_frame)} dwarfs")
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
                        wcs=self.wcs,
                        bbox=self.bbox,
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
            logger.info("not ingesting catalogs to change use self.run(ingest=True)")
            return self.injection_cats

    def generate_catalogs(self, multiproc=False):
        self.config["injection"]["mag_limit_band"]

        rows = self.dwarf_params_frame.to_dict(orient="records")
        r_len = len(rows)
        gen_args = list(
            zip(
                rows,
                r_len * [self.wcs],
                r_len * [self.bbox],
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

    def get_data_ids(self, tract=9615, patch=3):
        """
        Get data IDs based on band, skymap, tract, and patch information.

        Parameters
        ----------
        tract : int, optional
            Tract number. Default is 9615.
        patch : int, optional
            Patch number. Default is 3.
        """
        self.dataid_dict = {}
        for band in self.config["pipelines"]["bands"]:
            self.dataid_dict[band] = {
                "band": band,
                "skymap": "hsc_rings_v1",
                "tract": tract,
                "patch": patch,
            }

    def get_coadds(self):
        """
        Fetch coadd data based on the specified bands and data IDs.
        """
        self.coadd_dict = {}
        for band in self.config["pipelines"]["bands"]:
            n_dataid = len(self.dataid_dict[band])
            # currently this only grabs one coadd per band?
            # if  n_dataid > 1:
            #     msg = f'{n_dataid} calexp data ids in {band}-band only using first one'
            #     logger.warning(msg)
            self.coadd_dict[band] = self.butler.get(
                "deepCoadd_calexp", dataId=self.dataid_dict[band]
            )
        self.wcs = self.coadd_dict["g"].getWcs()
        self.bbox = self.coadd_dict["g"].getBBox()

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


if __name__ == "__main__":
    """
    Main execution script.
    Parses input arguments, reads the configuration file,
    and runs the dwarf injection catalog creation.
    """
    config_dict = default_config_dict.copy()
    parser = argparse.ArgumentParser(
        description="Description of your script/functionality."
    )
    parser.add_argument("filename", help="name of config.yaml file")
    args = parser.parse_args()

    with open(args.filename, "r") as file:
        user_config_dict = yaml.safe_load(file)
        for key, value in user_config_dict.items():
            config_dict[key] = value

    creator = CreateDwarfInjectionCatalog(config_dict)
    creator.run(ingest=False)
