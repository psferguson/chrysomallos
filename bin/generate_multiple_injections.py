import argparse

import numpy as np
from tqdm import tqdm

from chrysomallos.injection import (
    CreateDwarfInjectionCatalog,
    DwarfParamSampler,
    InjectionGenerator,
)
from chrysomallos.utils import Config
from chrysomallos.utils.annotations import get_anotation_box


def gal_size(r_scale, dwarf_params_frame, wcs, bbox):
    row = dwarf_params_frame.iloc[0].copy()
    row["r_scale"] = r_scale
    result = get_anotation_box(
        wcs=wcs, bbox=bbox, dwarf_params=row, scaling_factor=0.5, central_density=1
    )
    return (
        np.log10(r_scale),
        result["x_max"] - result["x_min"],
        result["y_max"] - result["y_min"],
    )


def run_single_config(args):
    tract, patch, sampling_seed, config_dict = args
    # when using multiproc we want different random seeds
    # for each process
    np.random.seed(sampling_seed)
    
    # read in the config file
    config = Config(config_dict)

    # set the tract and patch
    config["pipelines"]["tract"] = tract
    config["pipelines"]["patch"] = patch
    config["sampling"]["generation_id"] = np.random.randint(0, 1000000)
    # sample dwarfs from 1 to 3
    if isinstance(config["sampling"]["n_dwarfs"], list):
        assert len(config["sampling"]["n_dwarfs"]) == 2
        ndwarfs = np.random.randint(config["sampling"]["n_dwarfs"][0], 
                                    config["sampling"]["n_dwarfs"][1])
    
    elif isinstance(config["sampling"]["n_dwarfs"], int):
        ndwarfs = config["sampling"]["n_dwarfs"]
    else:
        raise ValueError(
            f"n_dwarfs must be either an integer or a list {config['sampling']['n_dwarfs']}"
            )
    print(f"{tract=} {patch=} {ndwarfs=}")
    config["sampling"]["n_dwarfs"] = ndwarfs
    # here is where we would want to customize the output collection
    # if config["stamp"]["version"] is not None:
    #     version = config["stamp"]["version"]
    #     config["sampling"]["output_directory"] = config["sampling"][
    #         "output_directory"
    #     ].format(version=version)
    #     config["stamp"]["directory"] = config["stamp"]["directory"].format(
    #         version=version
    #     )
    #     config["stamp"]["title_format"] = f"hsc_stamp_r_{version}"
    # config["sampling"][
    #     "output_file"
    # ] = f"round_{version}_tract_{tract}_patch_{patch}_{config['sampling']['generation_id']}.csv"
    config["sampling"]["random_seed_sampling"] = sampling_seed
    
    sampler = DwarfParamSampler(config)
    # logger.info(f"generating params for {config['sampling']['n_dwarfs']} dwarfs")
    dwarf_params_frame, coadd_dict = sampler.run(write=True)

    creator = CreateDwarfInjectionCatalog(config, dwarf_params_frame)
    catalogs, coadd_dict = creator.run(
        ingest=False,
        multiproc=False,
    )
    

    injection_generator = InjectionGenerator(
        config=config,
        dwarf_params_frame=dwarf_params_frame,
        dwarf_catalogs=catalogs,
        coadd_dict=coadd_dict,
    )
    injection_generator.run()
    

def run_configs(stamp_list, multiproc=False):
    gen_args = stamp_list  # [(tract, patch, ) for tract in tracts for patch in patches]

    if multiproc:
        from multiprocessing import Pool

        multiproc = int(multiproc)
        processes = multiproc if multiproc > 0 else None
        with Pool(processes, maxtasksperchild=1) as p:
            for _ in tqdm(
                p.imap_unordered(run_single_config, gen_args), total=len(gen_args)
            ):
                pass
    else:
        [run_single_config(args) for args in gen_args]


if __name__ == "__main__":
    """
    Main execution script.
    Parses input arguments, reads the configuration file,
    and runs the dwarf injection catalog creation to generate postage stamps.
    """
    parser = argparse.ArgumentParser(
        description="Description of your script/functionality."
    )
    parser.add_argument("--config", "-c", help="name of config.yaml file", default=None)
    parser.add_argument(
        "--multiproc",
        "-mp",
        default=False,
        help="use multiprocessing",
    )
    parser.add_argument(
        "--tracts",
        "-t",
        default=None,
        help="which tracts to use comma separated list",
    )
    parser.add_argument(
        "--n_inject",
        "-ni",
        default=1,
        type=int,
        help="how many times to inject dwarfs in each tract"
             "patch combination",
    )

    args = parser.parse_args()
    
    patches = np.arange(1) 
    if args.tracts is None:
        tracts = [9615, 9697, 9813]
    else:
        tracts = [int(tract) for tract in args.tracts.split(",")]

    injection_list = []  # [(tract, patch, sampling_seed, config_dict)]
    for tract in tracts:
        for patch in patches:
            injection_list.append(
                (tract, patch, np.random.randint(0, 1000000), args.config)
            )
    # could be optimized to check that all bands exist in repo before adding to 
    # injection_list

    n_iterations = args.n_inject
    
    print(f"{n_iterations=}")

    modified_stamp_list = injection_list
    if n_iterations > 1:
        for iteration in range(n_iterations - 1):
            modified_stamp_list = modified_stamp_list + [
                (tract, patch, np.random.randint(0, 1000000), config)
                for (tract, patch, __, config) in injection_list
            ]

    run_configs(modified_stamp_list, multiproc=args.multiproc)
