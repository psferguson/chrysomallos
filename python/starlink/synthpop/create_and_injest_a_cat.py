import argparse
import astropy.table as atable
import astropy.units as u
import numpy as np
from tqdm import tqdm
import yaml


from lsst.daf.butler import Butler, CollectionType
from lsst.daf.butler.registry import MissingCollectionError
import lsst.source.injection as si

from starlink.synthpop import adopt_a_cat, massage_the_cat
from starlink.utils.log import logger


default_config_dict={
    "repo" : "/repo/main",
    "collection" : 'HSC/runs/RC2/w_2023_32/DM-40356',
    "bands" : ["g","r","i"],
    "dwarfs": []
    }

class DwarfConfig():
    """
    Configuration class for a dwarf galaxy. 
    Sets various properties of the dwarf galaxy, including its location, 
    age, metallicity, and other parameters.
    """
    # DEFAULT PARAMETERS 
    # ID of the injection cat
    id = None
    # X-coordinate of the object's center [pixel]
    x_cen = None
    # Y-coordinate of the object's center [pixel]
    y_cen = None
    # Age of the object [Gyr]
    age = 10.0
    # Metallicity of the object [Fe/H]
    feh = -2.0
    # Mass of the object [Msolar]
    mass = 5.0e5
    # Distance to the object [Mpc]
    dist = 2.0 
    # Scale radius of the object [pc]
    r_scale = 100 
    # Ellipticity of the object (0-1)
    ellip = 0
    # Orientation angle of the object
    theta = 0
    # Sersic index of the object (shape parameter)
    n = 1
    # Magnitude limit to inject for the object
    mag_limit = 27
    # Band in which the magnitude limit is defined
    mag_limit_band = "LSST_g"
    # Seed for random number generation (for reproducibility)
    random_seed = None
    
    def __init__(self, dwarf_dict):
        """
        Initialize the class with the given dwarf dictionary.
        
        Parameters
        ----------
        dwarf_dict : dict
            Dictionary containing the properties and attributes of a dwarf.
        """
        # update new values
        for key, value in dwarf_dict.items():
            setattr(self, key, value)
            
        

class CreateDwarfInjectionCatalog():
    """
    Class to create an injection catalog for dwarf galaxies.
    Processes configuration data, fetches coadd from butler, 
    and creates catalogs ready for injection, then ingests into butler.
    """
    def __init__(self, config_dict=default_config_dict):
        """
        Initialize the class with the given configuration dictionary.
        
        Parameters
        ----------
        config_dict : dict, optional
            Configuration dictionary with parameters and settings.
            Default is `default_config_dict`.
        """
        self.config_dict = config_dict
        for key, value in self.config_dict.items():
            setattr(self, key, value)
    
    def run(self):
        """
        Main method to run the dwarf injection catalog creation process.
        """
        self.get_data_ids()
        self.butler = Butler(self.repo, collections=self.collection)
        self.get_coadds()
        self.dwarf_cats=[]
        # generate xy catalog for each dwarf
        for dwarf_config in self.dwarf_configs:
            self.dwarf_cats.append(
                adopt_a_cat(wcs = self.wcs, 
                            bbox = self.bbox, 
                            r_scale = dwarf_config.r_scale, 
                            n = dwarf_config.n, 
                            ellip = dwarf_config.ellip,
                            theta = dwarf_config.theta, 
                            random_seed = dwarf_config.random_seed,
                )
            )
       
        logger.info(f"generated catalogs for {len(self.dwarf_configs)} dwarfs")
        # convert to ssi format and concatenate to single catalog per band
        self.injection_cats = {}
        for band in self.bands:
            self.injection_cats[band] = []
        for i, cat in enumerate(self.dwarf_cats):
            for band in self.bands: 
                self.injection_cats[band].append(
                    massage_the_cat(cat_inp=cat, 
                                    injection_maglim=self.dwarf_configs[i].mag_limit,
                                    band_for_injection = band,
                                    xcen = self.dwarf_configs[i].x_cen, 
                                    ycen = self.dwarf_configs[i].y_cen, 
                                    wcs = self.wcs,
                                    bbox = self.bbox,
                                    r_scale=self.dwarf_configs[i].r_scale,
                                    dist = self.dwarf_configs[i].dist
                                    )
                )
        for band in self.bands: 
            self.injection_cats[band] = atable.vstack(self.injection_cats[band])
        
        self.ingest_injection_catalogs(si_input_collection = self.inject_cat_collection, 
                                         catalogs = self.injection_cats, 
                                         bands = self.bands
        )
        
        
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
        for band in self.bands:
            self.dataid_dict[band] = {'band': band, 'skymap': 'hsc_rings_v1', 'tract': tract, 'patch': patch}
    
    def get_coadds(self):
        """
        Fetch coadd data based on the specified bands and data IDs.
        """
        self.coadd_dict = {}
        for band in self.bands: 
            n_dataid=len(self.dataid_dict[band])
            # currently this only grabs one coadd per band? 
            # if  n_dataid > 1:
            #     msg = f'{n_dataid} calexp data ids in {band}-band only using first one'
            #     logger.warning(msg)
            self.coadd_dict[band] = self.butler.get('deepCoadd_calexp', dataId=self.dataid_dict[band])
        self.wcs = self.coadd_dict["g"].getWcs()
        self.bbox = self.coadd_dict['g'].getBBox()
    
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
        writeable_butler =  Butler(self.repo, writeable=True)
        try:
            writeable_butler.removeRuns([si_input_collection])
        except MissingCollectionError:
            logger.info("Writing into a new RUN collection")
            pass
        else:
            logger.info("Prior RUN collection located and successfully removed")
        _ = writeable_butler.registry.registerCollection(si_input_collection, type=CollectionType.RUN)
        
        self.injection_cat_refs = {}
        for band in bands:
            self.injection_cat_refs[band] = si.ingest_injection_catalog(
                writeable_butler=writeable_butler,
                table=catalogs[band],
                band=band,
                output_collection=si_input_collection,
                dataset_type_name='injection_catalog'
            )
    
if __name__ == "__main__":
    """
    Main execution script. 
    Parses input arguments, reads the configuration file, 
    and runs the dwarf injection catalog creation.
    """
    config_dict = default_config_dict.copy()
    parser = argparse.ArgumentParser(description="Description of your script/functionality.")
    parser.add_argument("filename", help="name of config.yaml file")
    args = parser.parse_args()
    
    with open(args.filename, "r") as file:
        user_config_dict = yaml.safe_load(file)
        for key, value in user_config_dict.items():
            config_dict[key] = value
        config_dict["dwarf_configs"]=[]
        for dwarf in config_dict["dwarfs"]:
            config_dict["dwarf_configs"].append(DwarfConfig(dwarf))
        

    creator=CreateDwarfInjectionCatalog(config_dict)
    creator.run()