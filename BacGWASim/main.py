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
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=52)
    parser = argparse.ArgumentParser(
        prog="BacGWASim", description="Description", 
        formatter_class=formatter
    )
    parser.add_argument(
        "--version", action="version",
        version="%(prog)s version {version}".format(version=__version__)
    )

    # Arg group - Genome simulation
    gen_group = parser.add_argument_group("Genome simulation parameters")
    gen_group.add_argument(
        "--num-species",
        default=None, metavar="INT", type=int,
        help="Number of samples in the simulated population"
    )
    gen_group.add_argument(
        "--genome-length",
        default=None, metavar="INT", type=int,
        help="Length of the genome (bp)"
    )
    gen_group.add_argument(
        "--mutation-rate",
        default=None, metavar="[0-1[", type=float,
        help="Mutation rate"
    )
    gen_group.add_argument(
        "--recomb-rate",
        default=None, metavar="[0-1[", type=float,
        help="Recombination rate"
    )
    gen_group.add_argument(
        "--maf",
        default=None, metavar="[0-1[", type=float,
        help="Minor allele frequency threshold of rare alleles to be discarded"
    )
    gen_group.add_argument(
        "--num-var",
        default=None, metavar="INT", type=int,
        help="Number of simulated variants. If '-1', variant number will be solely a function of mutation rate"
    )

    # Arg group - Phenotype simulation
    phen_group = parser.add_argument_group("Phenotype simulation parameters")
    phen_group.add_argument(
        "--phen-type",
        default=None, choices=["cc", "quant"],
        help="Type of simulated phenotype. 'cc':binary case-control, 'quant': quantitative"
    )
    phen_group.add_argument(
        "--num-causal-var",
        default=None, metavar="INT", type=int,
        help="Number of causal markers"
    )
    phen_group.add_argument(
        "--causal-maf-min",
        default=None, type=float,
        help="Minimum Minor Allele Frequency (MAF) of causal markers"
    )
    phen_group.add_argument(
        "--causal-maf-max",
        default=None, type=float,
        help="Maximum Minor Allele Frequency (MAF) of causal markers"
    )
    phen_group.add_argument(
        "--causal-ld-max",
        default=None, type=float,
        help="Maximum permitted R2 score between pairs of causal markers in window size of 1000 candidate causal markers meeting --causal-maf-min and --causal-maf-max thresholds"
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

    # Extracting information from the config file
    if args.config:
        print("Processing args from config file.")
        for key in args.config:
            if key in config:
                config[key] = args.config[key]

    # Using cli args (priority)
    for arg in [arg for arg in vars(args) if arg is not "config"]:
        if args[arg] is not None:
            config[arg] = args[arg]

    # Running snakemake
    snakemake.snakemake(
        snakefile=snakefile_path,
        config=config,
        forceall=True,
    )


if __name__ == '__main__':
    main()
