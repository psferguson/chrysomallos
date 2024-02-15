import argparse

from chrysomallos.injection import DwarfParamSampler
from chrysomallos.utils import Config, logger

if __name__ == "__main__":
    """
    Generates a dwarf parameters based on specified configurations.
    These can then be injected into coadd images.
    """
    parser = argparse.ArgumentParser(
        description="Generates a dwarf parameters based on specified configurations."
    )
    parser.add_argument(
        "--config",
        "-c",
        help="Path to the config.yaml file. Defaults to an example configuration if not specified.",
        default=None,
    )
    parser.add_argument(
        "--write",
        "-w",
        default=False,
        action="store_true",
        help="Enables writing the sampled dwarfs to a CSV or FITS file.",
    )
    args = parser.parse_args()
    if args.config is None:
        config_dict = "./python/chrysomallos/config/example_stamps.yaml"
    else:
        config_dict = args.config

    config = Config(config_dict)

    sampler = DwarfParamSampler(config)

    generated_parameters, coadd_dict = sampler.run(write=args.write)

    logger.info("Sampling Completed")
    print(generated_parameters)
