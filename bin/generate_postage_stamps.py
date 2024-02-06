import argparse
import os

from chrysomallos.injection import CreateDwarfInjectionCatalog, DwarfParamSampler
from chrysomallos.utils import Config, logger, read_data_file

if __name__ == "__main__":
    """
    Main execution script.
    Parses input arguments, reads the configuration file,
    and runs the dwarf injection catalog creation.
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
    if args.config is None:
        config_dict = "./python/chrysomallos/config/example_stamps.yaml"
    else:
        config_dict = args.config

    config = Config(config_dict)
    if os.path.exists(config["sampling"]["output_file"]) & ~force_sampling:
        logger.info("reading dwarf parameters from file")
        dwarf_params_frame = read_data_file(config["sampling"]["output_file"])
        logger.info(f"parameters for {len(dwarf_params_frame)} dwarfs")
    else:
        sampler = DwarfParamSampler(config)
        logger.info(f"generating params for {config['sampling']['n_dwarfs']} dwarfs")
        dwarf_params_frame = sampler.run()

    creator = CreateDwarfInjectionCatalog(config, dwarf_params_frame)
    catalogs = creator.run(
        ingest=False,
    )

    logger.info("Done")
