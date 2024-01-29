"""
Classes for storing and updating config for postage stamps and grid injections
"""
from dataclasses import dataclass, field
from inspect import signature
from typing import List

import numpy as np


@dataclass
class BaseConfig:
    # currently only support DEEP COADD, should we change for visit?
    pipelines: dict = field(
        default_factory=lambda: {
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
        }
    )

    # ranges to sample over
    sampling: dict = field(
        default_factory=lambda: {
            # sampler to use for parameter generation
            # will default to this unless specified
            "default_sampler": np.random.uniform,
            "n_dwarfs": np.nan,
            "m_v": [np.nan],
            "sb": [np.nan],
            "ellip_range": [np.nan],
            "x_cen": [np.nan],
            "y_cen": [np.nan],
        }
    )
    injection: dict = field(default_factory=lambda: {"mag_lim": 29.0})

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


@dataclass
class InjectConfig(BaseConfig):
    ingest: bool = True
    # generator: grid_params


@dataclass
class StampConfig(BaseConfig):
    ingest: bool = False
    n_dwarfs: int = 0
    dist_range: List[float] = field(default_factory=lambda: [2, 3])
    Mv_range: List[float] = field(default_factory=lambda: [-5, -11])
    sb_range: List[float] = field(default_factory=lambda: [24, 27])
    ellip_range: List[float] = field(default_factory=lambda: [0, 0.5])
    x_cen_range: List[float] = field(default_factory=lambda: [300, 3800])
    y_cen_range: List[float] = field(default_factory=lambda: [300, 3800])
    # generator: sample_params
    stamp_x_size: int = 600
    stamp_y_size: int = 600
