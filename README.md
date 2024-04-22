
# Stochastic Analysis of Infinite-State Probabilistic Models

The framework takes PRISM models of CRNs and a CSL property of the form $P=? [\textrm{true   } \textrm{U}<=\textrm{T   }  \textrm{X}=\theta]$ as input, and generates and iteratively increasing lower-bound for the probability of the property on the CRN model.
At each iteration, the lower-bound probabilty and the size of the state-space on which this probability is computed is printed to the console.

# Installation

- Prerequisites

	- Stormpy

		- Stormpy is a set of python bindings for the STORM’s C++ API. The framework uses this library for building models, and checking the probability of properties on the built models. Stormpy can be installed following the instructions showed [here](https://moves-rwth.github.io/stormpy/installation.html#installation-steps).

	- STORM

		- Stormpy installation requires STORM to be available on the system. The framework additionally uses the  `storm-conv` binary to convert PRISM models to JANI models. `storm-conv` binary will be built following the standard installation of STORM. The instructions to compile STORM can be found [here](https://www.stormchecker.org/documentation/obtain-storm/build.html).

	- z3-py

		- This framework uses z3 and its python bindings z3-py as the underlying SMT-solver. z3-py can added to the python interpreter by running the following command: `pip install z3-solver` . More details on installation of z3-py can be found [here](https://github.com/Z3Prover/z3?tab=readme-ov-file#z3-bindings).

- ==Note: Both Stormpy and z3-py should be installed in the same python virtual environment.==

- Cloning the python scripts
	
	- After prerequisite libraries are installed, the python scripts can be cloned and accessed by following commands:
		 ```bash
		git clone --depth 1 --branch qest2024 https://github.com/fluentverification/bmc_counterexample.git
		cd bmc_counterexample
		```

### Docker Image

Alternatively, a docker container with all the prerequisites libraries and case studies can be downloaded from the docker hub:

```bash

To  be  released  soon.

```

# Running

The framework accepts CTMC models written in [PRISM modeling language](https://www.prismmodelchecker.org/manual/ThePRISMLanguage/Introduction) as input. Currently, the framework only accepts a restricted class of models where:

1. All commands in the model are labelled.

2. There is no module-renaming in the model.

3. Each command in the model must have exactly one update.

In addition to the PRISM model, the JANI model of the input CTMC must be passed to the framework as well. The JANI model of the case-studies in this repository are available in the directory containing each CRN's PRISM model. The following command uses the `storm-conv` binary to convert the PRISM-Model.sm which is a PRISM model to a JANI model:
```bash
torm-conv --prism ./PRISM-Model.sm --prop "Property of Interest" --tojani JANI-Model.jani -pc
```
To run the framework, a JSON file containing the necessary parameters should be passed to the python interpreter. This JSON file should follow this format:

```json
{
"model_name" : (string)  arbitrary  name  given  to  the  model,
"model_path" : (string)  the  path  to  the  PRISM  model,
"jani_path" : (string)  the  path  to  the  jani  model,
"csl_property" : (string)  CSL  property  that  is  going  to  be  checked,
"target_variable" : (string)  the  variable  of  interest  in  the  CSL  property,
"target_value" : (string)  the  value  of  the  variable  of  interest  in  the  CSL  property
}
```

To run the framework, this JSON file should be simply passed as an argument to the [Main.py](http://Main.py) script:

```bash
python  Main.py  file.json
```  

# Case Studies  

### Enzymatic Futile Cycle

Enzymatic Futile Cycle is a chemical reaction network consisting of 6 chemical species reacting through 6 reaction channels.

Link to the PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/qest2024/CRNs/enzymatic_futile_cycle/enzym_unb.sm)

Link to the JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/qest2024/CRNs/enzymatic_futile_cycle/enzymatic_futile_cycle.json)

The following command runs the framework on this model:

```bash
python  Main.py  path-to-json.json
```

### Motility Regulation

Motility regulation is a chemical reaction network of a gene regulatory network consisting of 9 chemical species reacting through 12 reaction channels.

Link to PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/qest2024/CRNs/motility_regulation/motility_unb.sm)

Link to JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/qest2024/CRNs/motility_regulation/motility_regulation.json) 

The following command runs the framework on this model:

```bash
python  Main.py  path-to-json.json
```

### Yeast Polarization

Yeast Polarization is a chemical reaction network model consisting of 7 chemical species reacting through 8 reaction channels.

Link to PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/qest2024/CRNs/yeast_polarization/yeast_unb.sm) 

Link to JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/qest2024/CRNs/yeast_polarization/yeast_polarization.json)

The following command runs the framework on this model:

```bash
python  Main.py  path-to-json.json
```

### Genetic Circuit0x8E  

Genetic Circuit0x8E is a chemical reaction network of a genetic circuit model consisting of 18 chemical species reacting through 15 reaction channels.

Link to PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/qest2024/CRNs/circuit0x8E/Circuit0x8E_100to111_unb.sm)

Link to JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/qest2024/CRNs/circuit0x8E/circuit0x8E.json)

The following command runs the framework on this model:

```bash
python  Main.py  path-to-json.json
```

# Contact

Email: mahmadi@usf.edu

SEES lab website: [Link](https://sees-usf.github.io/)
