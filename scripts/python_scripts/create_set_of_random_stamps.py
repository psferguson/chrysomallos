from dataclasses import dataclass, field
from inspect import signature
from typing import List

import numpy as np
import pandas as pd
import yaml
from lsst.daf.butler import Butler

from chrysomallos.utils import mstar_from_absmag, sb_mv_to_rh

# given:
# - tract
# - config
# - nstamps inject/blank


# def run_generate_photometry(spec_files,download_directory,output_directory,log_directory,multiproc=False):
#     #logger = logging.getLogger()
#     if isinstance(spec_files, str):
#         spec_files = [spec_files]

#     N=len(spec_files)
#     #logging.debug("Generating Photometry for %s files..."%N)

#     inputs = zip(spec_files,N*[download_directory],N*[output_directory],N*[log_directory])


#     if multiproc:
#         from multiprocessing import Pool
#         processes = multiproc if multiproc > 0 else None
#         p = Pool(processes,maxtasksperchild=1)
#         out = p.map(generate_photometry,inputs)
#     else:
#         out = [generate_photometry(arg) for arg in inputs]


# Create base class for full catalog injection
@dataclass
class stampConfig:
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

    # stamp info
    stamp_x_size: int = 600
    stamp_y_size: int = 600

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


def sample_dwarf_params(stamp_config):
    distances = np.random.uniform(
        stamp_config.dist_range[0], stamp_config.dist_range[1], stamp_config.n_dwarfs
    )
    Mvs = np.random.uniform(
        stamp_config.Mv_range[0], stamp_config.Mv_range[1], stamp_config.n_dwarfs
    )
    sbs = np.random.uniform(
        stamp_config.sb_range[0], stamp_config.sb_range[1], stamp_config.n_dwarfs
    )
    r_scales = sb_mv_to_rh(sb=sbs, M_v=Mvs, distance=distances * 1e6)

    ellips = np.random.uniform(
        stamp_config.ellip_range[0], stamp_config.ellip_range[1], stamp_config.n_dwarfs
    )
    # assuming all cens tract coordinates
    x_cens = np.random.uniform(
        stamp_config.x_cen_range[0], stamp_config.x_cen_range[1], stamp_config.n_dwarfs
    )

    y_cens = np.random.uniform(
        stamp_config.y_cen_range[0], stamp_config.y_cen_range[1], stamp_config.n_dwarfs
    )

    # tract, patch, id, dist, x_cen, y_cen, r_scale, sb, mag_limit, m_v, mass, random_seed
    rng = np.random.default_rng()
    seeds = rng.integers(-2147483647, high=2147483646, size=stamp_config.n_dwarfs)
    param_frame = pd.DataFrame(
        {
            "tract": stamp_config.tract,
            "patch": stamp_config.patch,
            "id": np.arange(stamp_config.n_dwarfs),
            "dist": distances,
            "x_cen": x_cens,
            "y_cen": y_cens,
            "r_scale": r_scales,
            "sb": sbs,
            "Mv": Mvs,
            "mass": mstar_from_absmag(Mvs),
            "mag_limit": stamp_config.mag_lim,
            "random_seed": seeds,
        }
    )
    return param_frame


if __name__ == "__main__":
    import argparse

    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-sc", "--stamp_config", default=None)
    parser.add_argument("-if", "--injection_file", default=None)
    parser.add_argument("-t", "--tract", type=int, default=None)
    parser.add_argument("-p", "--patch", type=int, default=None)
    parser.add_argument("-n", "--n_dwarfs", type=int, default=None)

    args = parser.parse_args()

    # generation config
    if args.stamp_config is None:
        stamp_config = stampConfig()
    else:
        with open(args.stamp_config, "r") as file:
            stamp_config_dict = yaml.safe_load(file)
        stamp_config = stampConfig.from_kwargs(**stamp_config_dict)

    # update config based on argparse
    if args.tract is not None:
        stamp_config.tract = args.tract
    if args.patch is not None:
        stamp_config.patch = args.patch
    if args.n_dwarfs is not None:
        stamp_config.n_dwarfs = args.n_dwarfs

    stamp_config.n_dwarfs = 500
    # setup butler
    butler = Butler(stamp_config.repo, collections=stamp_config.collection)

    # get wcs
    image_dict = {}
    for band in stamp_config.bands:
        dataid = {
            "band": band,
            "skymap": stamp_config.skymap,
            "tract": stamp_config.tract,
            "patch": stamp_config.patch,
        }

        image_dict[band] = butler.get("deepCoadd_calexp", dataId=dataid)

    wcs = image_dict[stamp_config.bands[0]].getWcs()
    bbox = image_dict[stamp_config.bands[0]].getBBox()

    stamp_config.get_pix_cen_ranges(bbox)

    # generate dwarf properties and then save to csv
    dwarf_param_frame = sample_dwarf_params(stamp_config)
    # To Do: put runid in stamp config
    dwarf_param_frame.to_csv("./dwarf_params_v1.csv", index=False)


# generate all injection catalogs for that tract for each catalog we need
