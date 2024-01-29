"""
Classes for storing and updating config for postage stamps and grid injections
"""
from dataclasses import dataclass, field
from inspect import signature
from typing import List


@dataclass
class BaseConfig:
    # DEEP COADD
    # repository path
    repo: str = "/repo/main"
    # collection string
    collection: str = "HSC/runs/RC2/w_2023_32/DM-40356"
    # bands for injection
    bands: List[str] = field(default_factory=lambda: ["g", "r", "i"])
    # skymap
    skymap: str = "hsc_rings_v1"
    # tract to use
    tract: int = 9615
    # patch to use
    patch: int = 3

    # ranges to sample over
    n_dwarfs: int = 0
    dist_range: List[float] = field(default_factory=lambda: [2, 3])
    Mv_range: List[float] = field(default_factory=lambda: [-5, -11])
    sb_range: List[float] = field(default_factory=lambda: [24, 27])
    ellip_range: List[float] = field(default_factory=lambda: [0, 0.5])
    x_cen_range: List[float] = field(default_factory=lambda: [300, 3800])
    y_cen_range: List[float] = field(default_factory=lambda: [300, 3800])
    # ToDo: sersic index
    # injection info
    mag_lim: float = 29

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
    # generator: sample_params
    stamp_x_size: int = 600
    stamp_y_size: int = 600
