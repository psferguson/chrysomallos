import argparse
import astropy.table as atable
import astropy.units as u
import numpy as np
from tqdm import tqdm
import yaml
from dataclasses import dataclass


from lsst.daf.butler import Butler, CollectionType
from lsst.daf.butler.registry import MissingCollectionError
import lsst.source.injection as si

from starlink.synthpop import default_config_dict, CreateDwarfInjectionCatalog 
from starlink.utils.log import logger

if __name__ == "__main__":
    """
    Main execution script. 
    Parses input arguments, reads the configuration file, 
    and runs the dwarf injection catalog creation.
    """
    config_dict = default_config_dict.copy()
    parser = argparse.ArgumentParser(description="Description of your script/functionality.")
    parser.add_argument("filename", help="name of config.yaml file")
    parser.add_argument("--i","-ingest", default=False, action="store_true",
                        help="ingest catalogs into butler")
    
    args = parser.parse_args()
    
    with open(args.filename, "r") as file:
        user_config_dict = yaml.safe_load(file)
        for key, value in user_config_dict.items():
            config_dict[key] = value
        

    creator=CreateDwarfInjectionCatalog(config_dict)
    creator.run(ingest=True)