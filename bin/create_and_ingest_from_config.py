import argparse

import yaml

from chrysomallos.synthpop import CreateDwarfInjectionCatalog, default_config_dict

if __name__ == "__main__":
    """
    Main execution script.
    Parses input arguments, reads the configuration file,
    and runs the dwarf injection catalog creation.
    """
    config_dict = default_config_dict.copy()
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

    args = parser.parse_args()
    ingest_bool = args.i

    with open(args.filename, "r") as file:
        user_config_dict = yaml.safe_load(file)
        for key, value in user_config_dict.items():
            config_dict[key] = value

    creator = CreateDwarfInjectionCatalog(config_dict)
    creator.run(ingest=ingest_bool)
