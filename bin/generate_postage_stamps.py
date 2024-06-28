import argparse
import os

from chrysomallos.injection import (
    CreateDwarfInjectionCatalog,
    DwarfParamSampler,
    PostageStampGenerator,
)
from chrysomallos.utils import Config, logger, read_data_file

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
    if args.config is None:
        config_dict = "./python/chrysomallos/config/example_stamps.yaml"
    else:
        config_dict = args.config

    config = Config(config_dict)
    if config["stamp"]["version"] is not None:
        version = config["stamp"]["version"]
        config["sampling"]["output_directory"] = config["sampling"][
            "output_directory"
        ].format(version=version)
        config["stamp"]["directory"] = config["stamp"]["directory"].format(
            version=version
        )
        config["stamp"]["title_format"] = f"hsc_stamp_r_{version}"
    if os.path.exists(config["sampling"]["output_file"]) & ~force_sampling:
        logger.info("reading dwarf parameters from file")
        dwarf_params_frame = read_data_file(config["sampling"]["output_file"])
        logger.info(f"parameters for {len(dwarf_params_frame)} dwarfs")
        coadd_dict = None
    else:
        sampler = DwarfParamSampler(config)
        logger.info(f"generating params for {config['sampling']['n_dwarfs']} dwarfs")
        dwarf_params_frame, coadd_dict = sampler.run()

    creator = CreateDwarfInjectionCatalog(config, dwarf_params_frame, coadd_dict)
    catalogs, coadd_dict = creator.run(
        ingest=False,
        multiproc=5,
    )

    stampler = PostageStampGenerator(
        config=config,
        dwarf_params_frame=dwarf_params_frame,
        dwarf_catalogs=catalogs,
        coadd_dict=coadd_dict,
    )

    stampler.generate_empty_stamps(config["stamp"]["n_empty"])
    stampler.run()
    #
    logger.info("Done")
