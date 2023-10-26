# Counterexample Generation for Chemical Reaction Network Models using Bounded Model Checking

## Requirements
- [PRISM](https://www.prismmodelchecker.org/) binary
- [Storm](https://www.stormchecker.org/documentation/obtain-storm/build.html) binary
- [z3](https://github.com/Z3Prover/z3) solver
- [z3 python binding](https://github.com/Z3Prover/z3#python)
- [stormpy](https://moves-rwth.github.io/stormpy/installation.html)

## Running
`python3 Main.py __path_to_json__`


Model's path, necessary arguments, and the path for the prism/storm binaries are specified in a json file. To add a new file following elements must be set: 
- model : path to model
- starting_bound : the length of the shortest counterexample(witness) trace. If not known, set to 0.
- prism_binary : location of the prism's binary
- csl_property : a CSL property of the form "P=? [true U<=t X=Value]"
- property_variable : The variable(species) of interest, for the above CSL it would be X.
- property_value : The value of the variable we are checking for in the CSL property. For the above property it would be Value.
