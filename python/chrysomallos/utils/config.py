"""
Classes for storing and updating config for postage stamps and grid injections
"""
import copy
import os
import pprint
from inspect import signature

import numpy as np
import yaml

default_dict = {
    "pipelines": {
        "repo": None,
        "input_collection": None,
        "bands": None,
        "skymap": None,
        "tract": np.nan,
        "patch": np.nan,
    },
    "sampling": {
        "n_dwarfs": np.nan,
        "random_seed": np.nan,
        "params": {
            "distance": np.nan,
            "m_v": np.nan,
            "surface_brightness": np.nan,
            "ellipticity": np.nan,
            "x_cen": np.nan,
            "y_cen": np.nan,
            "age": np.nan,
            "feh": np.nan,
            "stellar_mass": np.nan,
            "r_scale": None,
            "sersic_index": np.nan,
            "ra": None,
            "dec": None,
        },
    },
    "injection": {
        "dwarfs": None,
        "mag_limit": np.nan,
        "mag_limit_band": None,
        "ingest": np.nan,
    },
}
default_stamp_dict = {
    "pipelines": {
        # repository path
        "repo": "/repo/main",
        # collection string
        "input_collection": "HSC/runs/RC2/w_2023_32/DM-40356",
        # bands for injection
        "bands": ["g", "r", "i"],
        # skymap
        "skymap": "hsc_rings_v1",
        # tract to use
        "tract": 9615,
        # patch to use
        "patch": 3,
    },
    "sampling": {
        "n_dwarfs": 0,
        "params": {
            "distance": [2, 3, "linear"],
            "m_v": [-5, -11, "linear"],
            "surface_brightness": [24, 27, "linear"],
            "ellipticity": [0, 0.5, "linear"],
            "x_cen": [300, 3800, "linear"],
            "y_cen": [300, 3800, "linear"],
            "age": 10,
            "feh": -2.0,
            "stellar_mass": 500000,
            "sersic_index": 1,
        },
    },
    "injection": {
        "mag_limit": 29,
        "mag_limit_band": "LSST_g",
        "ingest": False,
        "stamp_size": [600, 600],
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

    # currently only support DEEP COADD, should we change for visit?
    # pipelines: dict = field(
    #     default_factory=lambda: {
    #         # repository path
    #         "repo": "/repo/main",
    #         # collection string
    #         "input_collection": "HSC/runs/RC2/w_2023_32/DM-40356",
    #         # bands for injection
    #         "bands": ["g", "r", "i"],
    #         # skymap
    #         "skymap": "hsc_rings_v1",
    #         # tract to use
    #         "tract": 9615,
    #         # patch to use
    #         "patch": 3,
    #     }
    # )

    # # parameter sampling config
    # sampling: dict = field(
    #     default_factory=lambda: {
    #         # sampler to use for parameter generation
    #         # will default to this unless specified
    #         "sampler": None,
    #         "n_dwarfs": np.nan,
    #         "m_v": [np.nan],
    #         "sb": [np.nan],
    #         "ellip_range": [np.nan],
    #         "x_cen": [np.nan],
    #         "y_cen": [np.nan],
    #     }
    # )
    # injection: dict = field(
    #     default_factory=lambda: {
    #         "mag_lim": np.nan,
    #         "stamp_size": [np.nan, np.nan],
    #         "ingest": False
    #         "inject_cat_collection": None
    #     }
    # )

    # ToDo: sersic index
    # injection info

    @classmethod
    def from_kwargs(cls, **kwargs):
        # fetch the constructor's signature
        cls_fields = {field for field in signature(cls).parameters}
        # set named values and ignore parameters we dont care about
        native_args, new_args = {}, {}
        for name, val in kwargs.items():
            if name in cls_fields:
                native_args[name] = val
            else:
                new_args[name] = val
        # create class with native args
        ret = cls(**native_args)

        return ret

    def get_pix_cen_ranges(self, bbox):
        self.x_cen_range = [
            bbox.beginX + int(self.stamp_x_size / 2),
            bbox.endX - int(self.stamp_x_size / 2),
        ]
        self.y_cen_range = [
            bbox.beginY + int(self.stamp_y_size / 2),
            bbox.endY - int(self.stamp_y_size / 2),
        ]

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
        for key in self.keys():
            for sub_key in self[key].keys():
                if (key == "sampling") and (sub_key == "params"):
                    # for sampling params we dont want to overwrite full dict
                    # instead just want to replace the configured values
                    for param in config[key][sub_key].keys():
                        self[key][sub_key][param] = config[key][sub_key][param]
                else:
                    self[key][sub_key] = config[key][sub_key]


# @dataclass
# class InjectConfig(BaseConfig):
#     ingest: bool = True
#     # generator: grid_params


# @dataclass
# class StampConfig(BaseConfig):
#     sampling: dict = field(
#         default_factory=lambda: {
#             # sampler to use for parameter generation
#             # will default to this unless specified
#             "sampler": "uniform_sampler",
#             "n_dwarfs": 0,
#             # distance range in Mpc
#             "dist": [2,3]
#             "m_v": [-5,-11],
#             "sb": [24,27],
#             "ellip": [0,0.5],
#             "x_cen": [300,3800],
#             "y_cen": [300,3800],
#         }
#     )
#     "mag_lim": 29.0
#     ingest: bool = False
#     n_dwarfs: int = 0
#     dist_range: List[float] = field(default_factory=lambda: [2, 3])
#     Mv_range: List[float] = field(default_factory=lambda: [-5, -11])
#     sb_range: List[float] = field(default_factory=lambda: [24, 27])
#     ellip_range: List[float] = field(default_factory=lambda: [0, 0.5])
#     x_cen_range: List[float] = field(default_factory=lambda: [300, 3800])
#     y_cen_range: List[float] = field(default_factory=lambda: [300, 3800])
#     # generator: sample_params
#     stamp_x_size: int = 600
#     stamp_y_size: int = 600
