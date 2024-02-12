"""
Classes for storing and updating config for postage stamps and grid injections
"""
import copy
import os
import pprint

import numpy as np
import yaml

from chrysomallos.utils.log import logger

__all__ = [
    "default_dict",
    "Config",
]
default_dict = {
    "pipelines": {
        "repo": None,
        "input_collections": None,
        "bands": None,
        "skymap": None,
        "tract": np.nan,
        "patch": np.nan,
    },
    "sampling": {
        # sample or grid
        "type": None,
        "n_dwarfs": np.nan,
        "random_seed_sampling": np.nan,
        "params": {
            "distance": None,
            "m_v": None,
            "surface_brightness": None,
            "ellipticity": None,
            "theta": None,
            "x_cen": None,
            "y_cen": None,
            "age": None,
            "feh": None,
            "stellar_mass": None,
            "r_scale": None,
            "n": None,
            "ra": None,
            "dec": None,
            "random_seed_injection": None,
        },
        "output_file": None,
    },
    "injection": {
        "dwarfs": None,
        "mag_limit": np.nan,
        "mag_limit_band": None,
        "ingest": False,
        "output_collection": None,
        "type": None,  # either "grid" or "stamp"
    },
    "stamp": {
        "directory": "./stamps",
        "size": None,
        "title_format": "test_stamp_d_{}_mv_{}_sb_{}_r_{}.png",
        "Q": 10,
        "stretch": 0.5,
        "minimum": 0,
    },
}


class Config(dict):
    """
    Configuration Object
    """

    def __init__(self, config, default=default_dict):
        """
        initialize a configuration object from filename or dictonary,
        can also merge with default configuration"""
        self.update(self._load(default))
        self._update_params(self._load(config))

        # To Do: validation
        # self._validate()

    def write(self, filename):
        """
        Write a copy of this config object.

        Parameters:
        -----------
        outfile : output filename

        Returns:
        --------
        None
        """
        ext = os.path.splitext(filename)[1]
        writer = open(filename, "w")
        if ext == ".py":
            writer.write(pprint.pformat(self))
        elif ext == ".yaml":
            writer.write(yaml.dump(self))
        else:
            writer.close()
            raise Exception("Unrecognized config format: %s" % ext)
        writer.close()

    def _load(self, config):
        """Load this config from an existing config

        Parameters:
        -----------
        config : filename, config object, or dict to load

        Returns:
        --------
        params : configuration parameters
        """
        if isinstance(config, str):
            self.filename = config
            params = yaml.safe_load(open(config))
        elif isinstance(config, Config):
            # This is the copy constructor...
            self.filename = config.filename
            params = copy.deepcopy(config)
        elif isinstance(config, dict):
            params = copy.deepcopy(config)
        elif config is None:
            params = {}
        else:
            raise Exception("Unrecognized input")

        return params

    def _update_params(self, config):
        for key in config.keys():
            if key not in self.keys():
                logger.info(f"key '{key}' is not in the default config")
                continue
            for sub_key in config[key].keys():
                if sub_key not in self[key].keys():
                    logger.info(
                        f"key ['{key}']['{sub_key}'] is not in the default config"
                    )
                    continue
                if (key == "sampling") and (sub_key == "params"):
                    # for sampling params we dont want to overwrite full dict
                    # instead just want to replace the configured values
                    for param in config[key][sub_key].keys():
                        self[key][sub_key][param] = config[key][sub_key][param]
                else:
                    self[key][sub_key] = config[key][sub_key]
