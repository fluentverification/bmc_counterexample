# Counterexample Generation for Chemical Reaction Network Models using Bounded Model Checking

## Requirements
- [PRISM](https://www.prismmodelchecker.org/) binary
- [z3](https://github.com/Z3Prover/z3) solver
- [z3 python binding](https://github.com/Z3Prover/z3#python)

## Running
models are defined in models.py script. Currently there are 4 CRNs described in this script: *Single Species Production and Degradation*, *Enzymatic Futile Cycle*, *Yeast Polarization* and *Motility Regulation*. New models can be appended to the script.

Necessary arguments and location to the prism binary are specified in a json file in json/ subdirectory. To add a new file following elements must be set: 
- model : name of the model in *models.py* script
- starting_bound : the length of the shortest counterexample(witness) trace. If not known, set to 0.
- prism_binary : location of the prism's binary
- csl_property : a CSL property of the form "P=? [true U<=t X=Value]"
- property_variable : The variable(species) of interest, for the above CSL it would be X.
- property_value : The value of the variable we are checking for in the CSL property. For the above property it would be Value.
- probability_bound : The threshold for which a counterexample is to be generated. Must be some value between 0 and 1. 
- model_check_step : larger values would result in less calls to prism, therefore less overhead.
- construct_path_bound : Maximum length of added sub-traces using scaffolding method.
- base : set to initial value of property_variable in the model,
- step : set to (property_value-base) if divide and conquer is not intended.

To run the program, simply pass one of the json files to the *bmc_cex.py* script.


`python3 bmc_cex.py json/single_species.json` 


The resulting graph would be saved in PRISM explicit format in three files in the results/model_name folder: *model_name.tra*, *model_name.sta*, *model_name.lab*. Console log will be saved in the same subdirectory as *mode_name.results*.