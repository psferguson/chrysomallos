import argparse
import os
import numpy as np

from chrysomallos.injection import (
    CreateDwarfInjectionCatalog,
    DwarfParamSampler,
    PostageStampGenerator,
)
from chrysomallos.utils import Config, logger, read_data_file


def run_single_config(args):
    tract,patch, sampling_seed = args


    config_dict = "./config/stamps_home_net.yaml"
    config = Config(config_dict)
    
    config["pipelines"]["tract"] = tract
    config["pipelines"]["patch"] = patch
    config["sampling"]["output_file"] = f"round_3_tract_{tract}_patch_{patch}.csv"
    config["sampling"]["random_seed_sampling"] = sampling_seed
    
    sampler = DwarfParamSampler(config)
    logger.info(f"generating params for {config['sampling']['n_dwarfs']} dwarfs")
    dwarf_params_frame, coadd_dict = sampler.run(write=True)

    creator = CreateDwarfInjectionCatalog(config, dwarf_params_frame)
    catalogs, coadd_dict = creator.run(
        ingest=False,
        multiproc=False,
    )

    postage_stamp_generator = PostageStampGenerator(
        config=config,
        dwarf_params_frame=dwarf_params_frame,
        dwarf_catalogs=catalogs,
        coadd_dict=coadd_dict,
    )

    postage_stamp_generator.run()
    postage_stamp_generator.generate_empty_stamps(config["stamp"]["n_empty"])
    

def run_configs(stamp_list, multiproc=False):

    gen_args = stamp_list#[(tract, patch, ) for tract in tracts for patch in patches]
    
    if multiproc:
        from multiprocessing import Pool

        processes = multiproc if multiproc > 0 else None
        p = Pool(processes, maxtasksperchild=1)
        p.map(run_single_config, gen_args)
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
        "--force_sampling",
        "-fs",
        default=False,
        action="store_true",
        help="write csv of sampled dwarfs",
    )
    args = parser.parse_args()
    force_sampling = args.force_sampling
    patches = [1,  28, 31, 34, 36, 42, 63, 64, 65, 67,73,79]
    tracts = [9615, 9697, 9813]
    from glob import glob

    flist = glob("deepCoadd_repo/HSC/runs/RC2/w_2023_32/DM-40356/20230819T003257Z/deepCoadd_calexp/*/*/*/*")
    directory = "deepCoadd_repo/HSC/runs/RC2/w_2023_32/DM-40356/20230819T003257Z/deepCoadd_calexp/"
    stamp_list = []
    for tract in tracts:
        for patch in np.arange(81):
            band_count = 0
            for band in ['g', 'r', 'i']:
                if os.path.exists(directory + "/" + str(tract) + "/" + str(patch) + "/" + band):
                    band_count += 1
            if band_count == 3:
                stamp_list.append((tract, patch, np.random.randint(0, 1000000)))

    run_configs(stamp_list, multiproc=36)
