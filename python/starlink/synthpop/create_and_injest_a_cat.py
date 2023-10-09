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
default_dwarf_dict = {
    "id": 0,
    "x_cen": None,
    "y_cen": None,
    "age": 10.0,
    "feh": -2.0,
    "mass": 5.0e6,
    "dist": 2.0,
    "r_scale": 100,
    "ellip": 0,
    "theta": 0,
    "n": 1,
    "mag_limit": 27,
    "mag_limit_band": "LSST_g",
    "random_seed": None
}
class DwarfConfig():
    def __init__(self, dwarf_dict):
        #set defaults:
        config_dict = default_dwarf_dict.copy()
        
        # update new values
        for key, value in dwarf_dict.items():
            config_dict[key] = value
       
        for key, value in config_dict.items():
            setattr(self, key, value)
        

class CreateDwarfInjectionCatalog():
    """ some stuff here """
    def __init__(self, config_dict=default_config_dict):
        self.config_dict = config_dict
        for key, value in self.config_dict.items():
            setattr(self, key, value)
    
    def run(self):
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
        self.dataid_dict = {}
        for band in self.bands:
            self.dataid_dict[band] = {'band': band, 'skymap': 'hsc_rings_v1', 'tract': tract, 'patch': patch}
    
    def get_coadds(self):
        self.coadd_dict = {}
        for band in self.bands: 
            n_dataid=len(self.dataid_dict[band])
            # if  n_dataid > 1:
            #     msg = f'{n_dataid} calexp data ids in {band}-band only using first one'
            #     logger.warning(msg)
            self.coadd_dict[band] = self.butler.get('deepCoadd_calexp', dataId=self.dataid_dict[band])
        self.wcs = self.coadd_dict["g"].getWcs()
        self.bbox = self.coadd_dict['g'].getBBox()
    
    def ingest_injection_catalogs(self, si_input_collection, catalogs, bands):
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