import os

def parsing_config(config):

    # Must be percentages
    are_percentages_exclusive = [
        "mutation_rate", "recomb_rate", "maf",
        "causal_maf_min", "causal_maf_max",
        "disease_prevalence", "ld_maf"
    ]

    are_percentages_inclusive = [
        "causal_ld_max", "heritability"
    ]

    for param in are_percentages_exclusive:
        is_percentage_exclusive(param, config[param])

    for param in are_percentages_inclusive:
        is_percentage_inclusive(param, config[param])

    # effect_size_odr is multiple of num_causal_var
    if config["num_causal_var"] % len(config["effect_size_odr"].split(',')):
        raise ValueError("--effect-size-odr list length must be a multiple of --num-causal-var")

    # num_var must be non-null int or -1
    if (config["num_var"] == 0) or (config["num_var"] < -1):
        raise ValueError("--number-var must be a non-null int, or -1")

    # case + control < num_species
    if config["phen_type"] == "cc":
        if (config["case"] + config["control"] > config["num_species"]):
           raise ValueError("--num-species must be bigger than the sum of --case and --control")

    # causal_maf_min must be smaller than causal_maf_max
    if config["causal_maf_min"] >= config["causal_maf_max"]:
        raise ValueError("--causal-maf-min must be smaller than --causal-maf-max")

    # plot_ld must be a bool
    if config["plot_ld"] not in [True, False]:
        raise ValueError("--plot-ld must be True or False")


def is_percentage_exclusive(param, val):
    if (val <= 0) or (val >= 1):
        raise ValueError(param + " value must be between 0 and 1 (exclusive, ]0-1[)")

def is_percentage_inclusive(param, val):
    if (val < 0) or (val > 1):
        raise ValueError(param + " value must be between 0 and 1 (inclusive, [0-1])")


def config2file(config, template):
    with open(template, "r") as file:
        template_data = file.read()

    for key, val in config.items():
        template_data = template_data.replace("$" + key.upper(), str(val))

    template_output = os.path.join(config["output_dir"], "configfile.yaml")
    with open(template_output, "w") as file:
        file.write(template_data)
