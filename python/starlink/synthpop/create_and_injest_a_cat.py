import numpy as np
import astropy.units as u

from lsst.daf.butler import Butler, DimensionUniverse, DatasetType, CollectionType
from lsst.daf.butler.registry import MissingCollectionError
import lsst.source.injection as si
from tqdm import tqdm
from starlink.synthpop import adopt_a_cat, massage_the_cat
from starlink.utils.log import logger


default_config_dict={
    "repo" : "/repo/main",
    "collection" : 'HSC/runs/RC2/w_2023_32/DM-40356',
    "bands" : ["g","r","i"]
    }


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
        logger.info(self.bbox)
        
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
        
if __name__ == "__main__":
    creator=CreateDwarfInjectionCatalog()
    creator.run()