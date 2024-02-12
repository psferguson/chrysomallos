import argparse

from chrysomallos.injection import DwarfParamSampler
from chrysomallos.utils import Config, logger

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
        "--write",
        "-w",
        default=False,
        action="store_true",
        help="write csv of sampled dwarfs",
    )
    args = parser.parse_args()
    if args.config is None:
        config_dict = "./python/chrysomallos/config/example_stamps.yaml"
    else:
        config_dict = args.config

    config = Config(config_dict)

    sampler = DwarfParamSampler(config)

    generated_parameters = sampler.run(write=args.write)

    logger.info("Done")
    print(generated_parameters)
