import astropy.table as atable
import lsst.source.injection as si
import numpy as np
from lsst.daf.butler import Butler, CollectionType
from lsst.daf.butler.registry import MissingCollectionError

from chrysomallos.injection import adopt_a_cat, massage_the_cat
from chrysomallos.utils import get_coadd_dict, logger


class CreateDwarfInjectionCatalog:
    """
    Creates a catalog for injecting dwarf galaxies into images.
    Reads configuration, fetches coadd data, and prepares catalogs
    for injection. Optionally, catalogs can be ingested into a butler
    repository.
    """

    def __init__(self, config, dwarf_params_frame, coadd_dict=None):
        """
        Initializes the catalog creator with configuration and dwarf
        parameters.

        Parameters
        ----------
        config : dict
            Configuration parameters for the catalog creation.
        dwarf_params_frame : DataFrame
            Parameters for each dwarf galaxy to be injected.
        coadd_dict : dict, optional
            Coadd information for the image areas of interest.
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
        Executes the catalog creation and ingestion process.

        if self.config['sampling']['type']=='grid':
            all generated catalogs are concatenated at the end
            to a dictionary with keys being the bands
        if self.config['sampling']['type']=='stamp':
            a dictonary will be returned with keys being the bands
            and the values will be a list of catalogs one for each dwarf

        Parameters
        ----------
        ingest : bool, default=False
            If True, ingests the generated catalogs into the butler.
        multiproc : bool, default=False
            Enables multiprocessing for catalog generation.

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
                            dist=self.dwarf_params_frame["dist"][i],
                            theta=self.dwarf_params_frame["theta"][i],
                            ellip=self.dwarf_params_frame["ellip"][i],
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
                        dist=self.dwarf_params_frame["dist"][i],
                        theta=self.dwarf_params_frame["theta"][i],
                        ellip=self.dwarf_params_frame["ellip"][i],
                    )
        else:
            raise Exception(
                f"unknown injection type: {self.config['injection']['type']}"
            )

        return self.injection_cats, self.coadd_dict

    def generate_catalogs(self, multiproc=False):
        """
        Generates individual catalogs for each dwarf galaxy.

        Parameters
        ----------
        multiproc : bool, default=False
            If True, uses multiprocessing to speed up the process.
        """
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
        """
        Generates a single catalog for a dwarf galaxy.

        Parameters
        ----------
        args : tuple
            Contains the parameters for the dwarf galaxy.
        """
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
            dist=row["dist"],
            r_scale=row["r_scale"],
            ellip=row["ellip"],
            theta=row["theta"],
            n=row["n_sersic"],
            m_v=row["m_v"],
            mag_limit=36,  # notice this is set to be deeper
            mag_limit_band=mag_limit_band,
            random_seed=random_seed,
        )
        return catalog

    def ingest_injection_catalogs(self, si_input_collection, catalogs, bands):
        """
        Ingests generated catalogs into a butler collection.

        Parameters
        ----------
        si_input_collection : str
            The collection name for ingestion.
        catalogs : dict
            The catalogs to ingest, keyed by band.
        bands : list
            The bands to ingest catalogs for.
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
