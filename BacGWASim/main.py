import argparse
import os
import snakemake
import sys
import yaml

from . import _program, __version__
from .arg_parsing import parsing_config, config2file


def is_valid_config_file(parser, arg):
    try:
        handle = open(arg, "r")
    except IOError:
        parser.error("File %s not found." % arg)
    
    file = open(arg)
    return yaml.load(file, Loader=yaml.FullLoader)


def str2bool(arg):
    if arg.lower() in ["true", "yes", "y", "1"]:
        return True
    elif arg.lower() in ["false", "no", "n", "0"]:
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



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
    parser.add_argument(
        "--config", dest="config",
        metavar="FILE", type=lambda x: is_valid_config_file(parser, x),
        help="Path to a config file", 
    )
    parser.add_argument(
        "--output-dir",
        metavar="DIR",
        help="Path to the output directory", 
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
        default=None, metavar="]0-1[", type=float,
        help="Mutation rate"
    )
    gen_group.add_argument(
        "--recomb-rate",
        default=None, metavar="]0-1[", type=float,
        help="Recombination rate"
    )
    gen_group.add_argument(
        "--maf",
        default=None, metavar="]0-1[", type=float,
        help="Minor allele frequency threshold of rare alleles to be discarded"
    )
    gen_group.add_argument(
        "--num-var",
        default=None, metavar="INT", type=int,
        help="Number of simulated variants. If '-1', variant number will be solely a function of mutation rate"
    )
    gen_group.add_argument(
        "--random-seed",
        default=None, metavar="INT", type=int,
        help="Random seed for reproducibility of results"
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
        default=None, metavar="]0-1[", type=float,
        help="Minimum Minor Allele Frequency (MAF) of causal markers"
    )
    phen_group.add_argument(
        "--causal-maf-max",
        default=None, metavar="]0-1[", type=float,
        help="Maximum Minor Allele Frequency (MAF) of causal markers"
    )
    phen_group.add_argument(
        "--causal-ld-max",
        default=None, metavar="[0-1]", type=float,
        help="Maximum permitted R2 score between pairs of causal markers in window size of 1000 candidate causal markers meeting --causal-maf-min and --causal-maf-max thresholds"
    )
    phen_group.add_argument(
        "--effect-size-odr",
        default=None, metavar="LIST",
        help="Effect sizes of causal markers (.odds ratios) (comma separated, must be a multiple of --num-causal-var)"
    )
    phen_group.add_argument(
        "--phen-replication",
        default=None, metavar="INT", type=int,
        help="Number of phenotype replication sets"
    )
    phen_group.add_argument(
        "--heritability",
        default=None, metavar="[0-1]", type=float,
        help="Heritability of phenotype"
    )
    phen_group.add_argument(
        "--disease-prevalence",
        default=None, metavar="]0-1[", type=float,
        help="Prevalence of phenotype"
    )
    phen_group.add_argument(
        "--case",
        default=None, type=float,
        help="In case of case-control binary phenotype simulation, number of case and control samples must be defined by 'case' and 'control' parameters"
    )
    phen_group.add_argument(
        "--control",
        default=None, type=float,
        help="In case of case-control binary phenotype simulation, number of case and control samples must be defined by 'case' and 'control' parameters"
    )
    phen_group.add_argument(
        "--simvis-dpi",
        default=None, type=int,
        help="Set the DPI for the simVis plot"
    )

    # Arg group - Linkage disequilibrium plotting
    plot_group = parser.add_argument_group("Linkage Disequilibrium plotting")
    plot_group.add_argument(
        "--plot-ld",
        default=None, metavar="BOOL", type=str2bool,
        help="Generate the LD plot"
    )
    plot_group.add_argument(
        "--snp-limit",
        default=None, metavar="INT", type=int,
        help="Number of SNPs randomly selected for plotting linkage map (Increasing this number will significantly increase computation time)"
    )
    plot_group.add_argument(
        "--ld-maf",
        default=None, metavar="]0-1[", type=float,
        help="Minimum Minor Allele Frequency of markers for LD plotting (Lower this values, it is more difficult to estimate accurate r2 values between pairs of markers leading to more noisy plot)"
    )

    # Arg group - Runetime parameters
    run_group = parser.add_argument_group("Runtime parameters")
    run_group.add_argument(
        "--cores",
        default=None, metavar="INT", type=int,
        help="Number of cores available for computations"
    )
    run_group.add_argument(
        "--latency-wait",
        default=3, metavar="INT", type=int,
        help="Time to wait (in sec) after a job to ensure that all files are present"
    )

    # Parsing args
    args = parser.parse_args()

    # Loading default config values
    with open(config_path) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    # Extracting information from the config file
    if args.config:
        for key in args.config:
            if key in config:
                config[key] = args.config[key]

    # Using cli args (priority)
    args = vars(args)
    for arg in [arg for arg in args if arg != "config"]:
        if args[arg] != None:
            config[arg] = args[arg]

    # Parsing the configuration
    parsing_config(config)

    # Saving the modified config file in the results
    config2file(config)

    # Running snakemake
    snakemake.snakemake(
        snakefile=snakefile_path,
        config=config,
        forceall=True,
        cores=config["cores"],
        latency_wait=config["latency_wait"],
    )


if __name__ == '__main__':
    main()
