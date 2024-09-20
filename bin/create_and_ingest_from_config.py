import argparse
import os
import yaml

from chrysomallos import CreateDwarfInjectionCatalog, default_dict, DwarfParamSampler, logger, read_data_file
from chrysomallos import get_coadd_dict

if __name__ == "__main__":
    """
    Main execution script.
    Parses input arguments, reads the configuration file,
    and runs the dwarf injection catalog creation.
    """
    config = default_dict.copy()
    parser = argparse.ArgumentParser(
        description="Description of your script/functionality."
    )
    parser.add_argument("filename", help="name of config.yaml file")
    parser.add_argument(
        "-ingest",
        "--i",
        default=False,
        action="store_true",
        help="ingest catalogs into butler",
    )
    parser.add_argument(
        "--force_sampling",
        "-fs",
        default=False,
        action="store_true",
        help="write csv of sampled dwarfs",
    )

    args = parser.parse_args()
    ingest_bool = args.i
    force_sampling = args.force_sampling

    with open(args.filename, "r") as file:
        user_config_dict = yaml.safe_load(file)
        for key, value in user_config_dict.items():
            config[key] = value

    # creator = CreateDwarfInjectionCatalog(config_dict)
    # creator.run(ingest=ingest_bool)

    # print('Running get_coadd_dict')
    # get_coadd_dict(None, config)
    # print('Finished running get_coadd_dict')
    if os.path.exists(config["sampling"]["output_file"]) & ~force_sampling:
        logger.info("reading dwarf parameters from file")
        dwarf_params_frame = read_data_file(config["sampling"]["output_file"])
        logger.info(f"parameters for {len(dwarf_params_frame)} dwarfs")
        coadd_dict = None
    else:
        sampler = DwarfParamSampler(config)
        # print('Jesus was way cool')
        logger.info(f"generating params for {config['sampling']['n_dwarfs']} dwarfs")
        dwarf_params_frame, coadd_dict = sampler.run()

    creator = CreateDwarfInjectionCatalog(config, dwarf_params_frame, coadd_dict)
    catalogs, coadd_dict = creator.run(
        ingest=config['injection']['ingest'],
    )
#        multiproc=5,
    # import pdb; pdb.set_trace()
