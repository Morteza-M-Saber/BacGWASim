import argparse
import os
import snakemake
import sys
import yaml

from . import _program, __version__


def is_valid_config_file(parser, arg):
    try:
        handle = open(arg, "r")
    except IOError:
        parser.error("File %s not found." % arg)
    
    file = open(arg)
    return yaml.load(file, Loader=yaml.FullLoader)


def main(sysargs=sys.argv[1:]):
    # Path
    bac_path = os.path.dirname(os.path.realpath(__file__))
    snakefile_path = os.path.join(bac_path, "Snakefile")
    config_path = os.path.join(bac_path, "configfile.yaml")

    # Argparse parser
    parser = argparse.ArgumentParser(prog="BacGWASim", description="Description")
    parser.add_argument(
        "--version", action="version",
        version="%(prog)s version {version}".format(version=__version__)
    )

    # Config yaml file as input
    parser.add_argument(
        "--config", dest="config", help="This is help", metavar="FILE",
        type=lambda x: is_valid_config_file(parser, x)
    )

    # Parsing args
    args = parser.parse_args()

    # Loading default config values
    with open(config_path) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    # # Extracting information from the config file
    # if args.config:
    #     print("Processing args from config file.")
    #     for key in args.config:
    #         if key in config:
    #             config[key] = args.config[key]

    # print(config)

    # Path
    # bacg_dir = os.path.abspath(os.path.dirname(__file__))

    print('Starting snakemake')

    # # Config
    # config = dict()

    # config["simbac_path"] =  "../dependencies/sim/SimBac/SimBac"
    # config["haploview"] =  "../dependencies/ld/Haploview.jar"
    # config["gcta"] =  "../dependencies/gcta/gcta_1.93.2beta/gcta64"

    snakemake.snakemake(
        snakefile=snakefile_path,
        config=config,
        forceall=True,
    )


if __name__ == '__main__':
    main()
