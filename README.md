# Stochastic Analysis of Infinite-State Probabilistic Models

## [Installation](https://github.com/fluentverification/bmc_counterexample/tree/qest2024?tab=readme-ov-file#Installation) <pre>    </pre> [Running](https://github.com/fluentverification/bmc_counterexample/tree/qest2024?tab=readme-ov-file#Running)<pre>    </pre>[Case Studies](https://github.com/fluentverification/bmc_counterexample/tree/qest2024?tab=readme-ov-file#case-studies)<pre>    </pre>[Contact](https://github.com/fluentverification/bmc_counterexample/tree/qest2024?tab=readme-ov-file#Contact)

# Installation

- Prerequisites
    - Stormpy
        - Stormpy is a set of python bindings for the STORM’s C++ API. The framework uses this library for building models, and checking the probability of properties on the built models.  Stormpy can be installed following the instructions showed [here](https://moves-rwth.github.io/stormpy/installation.html#installation-steps).
    - STORM
        - The framework uses `storm-conv` binary to convert PRISM models to JANI models. STORM as the back-end probabilistic model checker. The instructions to compile STORM can be found [here](https://www.stormchecker.org/documentation/obtain-storm/build.html). The `storm-conv` binary is   After installation, the path to `storm-conv` binary should be declared in the JSON file passed to the framework.
    
    - z3-py
        - This framework uses z3 and its python bindings z3-py as the underlying SMT-solver. z3-py can added to the python interpreter by running the following command:    `pip install z3-solver` . More details on installation of z3-py can be found [here](https://github.com/Z3Prover/z3?tab=readme-ov-file#z3-bindings).
    - Note: Both Stormpy and z3-py should be installed in the same python virtual environment.
- After installing the prerequisites, the python3 with the installed libraries can be used to run the framework.
    
    ```bash
	git clone --depth 1 --branch qest2024 https://github.com/fluentverification/bmc_counterexample.git
	cd bmc_counterexample
    
    ```
    

### Docker Image

Alternatively, a docker container with all the prerequisites libraries and case studies can be downloaded from the docker hub:

```bash
To be released soon.
```


# Running

The framework accepts CTMC models written in [PRISM modeling language](https://www.prismmodelchecker.org/manual/ThePRISMLanguage/Introduction) as input. Currently, the framework only accepts a restricted class of models where:

1. All commands in the model are labelled.
2. There is no module-renaming in the model.
3. Each command in the model must have exactly one update.

To run the framework, a JSON file containing the necessary parameters should be passed to the python interpreter. This JSON file should follow this format:

```json
{
	"model_name" : (string) arbitrary name given to the model,
	"model_path" : (string) the path to the PRISM model,
	"storm-conv" : (string) the path to the storm-conv binary,
	"csl_property_lb" : (string) CSL property that is going to be checked,
	"target_variable" : (string) the variable of interest in the CSL property,
	"target_value" : (string) the value of the variable of interest in the CSL property
}
```

To run the framework, this JSON file should be simply passed as an argument to the [Main.py](http://Main.py) script:

```bash
python Main.py file.json
```

# Case Studies

### Enzymatic Futile Cycle

Enzymatic Futile Cycle is a chemical reaction network that is defined by the following set of reactions:

$$
\begin{array}{lll}
    R_1 : \ \textrm{S}_1 + \textrm{S}_2 \xrightarrow{1.0} \textrm{S}_3,~~~ &
    R_2 : \ \textrm{S}_3 \xrightarrow{1.0} \textrm{S}_1 + \textrm{S}_2,~~~ \\
    R_3 : \ \textrm{S}_3 \xrightarrow{0.1} \textrm{S}_1 + \textrm{S}_5, &
    R_4 : \ \textrm{S}_4 + \textrm{S}_5 \xrightarrow{1.0} \textrm{S}_6,~~~ \\
    R_5 : \ \textrm{S}_6 \xrightarrow{1.0} \textrm{S}_4 + \textrm{S}_5,~~~ &
    R_6 : \ \textrm{S}_6 \xrightarrow{0.1} \textrm{S}_4 + \textrm{S}_2
\end{array}
$$

Where the initial population of species $(S_1, S_2, S_3, S_4, S_5, S_6)$ are:

$$
\textbf{x0} = [1, 50, 0, 1, 50, 0].
$$

Link to the PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/main/CAV/CRNs/enzymatic_futile_cycle/enzym_unb.sm)

Link to the JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/main/CAV/CRNs/enzymatic_futile_cycle/enzymatic_futile_cycle.json)

### Motility Regulation

Motility regulation is a chemical reaction network of a gene regulatory network consisting of 9 chemical species reacting through 12 reaction channels:

$$
\begin{array}{ll}
    R_1 : \ \textrm{codY} \xrightarrow{0.1} \textrm{codY} + \textrm{CodY}, \\
    R_2 : \ \textrm{CodY} \xrightarrow{0.0002} \emptyset, \\
    R_3 : \ \textrm{flache} \xrightarrow{1} \textrm{flache} + \textrm{SigD}, \\
    R_4 : \ \textrm{SigD} \xrightarrow{0.0002} \emptyset, \\
    R_5 : \ \textrm{SigD\_hag} \xrightarrow{1} \textrm{SigD} + \textrm{hag} + \textrm{Hag}, \\
    R_6 : \ \textrm{Hag} \xrightarrow{0.0002} \emptyset, \\
    R_7 : \ \textrm{SigD} + \textrm{hag} \xrightarrow{0.01} \textrm{SigD\_hag}, \\
    R_8 : \ \textrm{SigD\_hag} \xrightarrow{0.1} \textrm{SigD} + \textrm{hag}, \\
    R_9 : \ \textrm{CodY} + \textrm{flache} \xrightarrow{0.02} \textrm{CodY\_flache},~~ \\
    R_{10} : \ \textrm{CodY\_flache} \xrightarrow{0.1} \textrm{CodY} + \textrm{flache}, \\
    R_{11} : \ \textrm{CodY} + \textrm{hag} \xrightarrow{0.01} \textrm{CodY\_hag}, \\
    R_{12} : \ \textrm{CodY\_hag} \xrightarrow{0.1} \textrm{CodY} + \textrm{hag} \\
\end{array}
$$

where the initial populations of the species $(\textrm{codY}, \textrm{CodY}, \textrm{flache}, \textrm{SigD}, \textrm{SigD\_hag}, \textrm{hag}, \textrm{Hag}, \textrm{CodY\_flache}, \textrm{CodY\_hag})$

are:

$$
\bm{x_0} = [1, 10, 1, 10, 1, 1, 10, 1, 1].
$$

Link to PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/main/CAV/CRNs/motility_regulation/motility_unb.sm)

Link to JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/main/CAV/CRNs/motility_regulation/motility_regulation.json)

### Yeast Polarization

Yeast Polarization is a chemical reaction network model consisting of 7 chemical species reacting through 8 reaction channels:

$$
\begin{array}{lll}
    R_1 : \ \emptyset \xrightarrow{0.0038} \textrm{R}, &
    R_2 : \ \textrm{R} \xrightarrow{4.00\times 10^{-4}} \emptyset,~~~ \\
    R_3 : \ \textrm{L} + \textrm{R} \xrightarrow{0.042} \textrm{RL} + \textrm{L}, &
    R_4 : \ \textrm{RL} \xrightarrow{0.0100} \textrm{R},~~~ \\
    R_5 : \ \textrm{RL} + \textrm{G} \xrightarrow{0.011} \textrm{G}_\textrm{a} + \textrm{G}_{\textrm{bg}},~~~ &
    R_6 : \ \textrm{G}_\textrm{a} \xrightarrow{0.100} \textrm{G}_\textrm{d}, \\
    R_7 : \ \textrm{G}_\textrm{d} + \textrm{G}_{\textrm{bg}} \xrightarrow{1.05\times 10^{3}} \textrm{G},~~~ &
    R_8 : \ \emptyset \xrightarrow{3.21} \textrm{RL} & \\
\end{array}
$$

Link to PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/main/CAV/CRNs/yeast_polarization/yeast_unb.sm)

Link to JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/main/CAV/CRNs/yeast_polarization/yeast_polarization.json)

### Genetic Circuit0x8E


# Contact

Email: mahmadi@usf.edu

SEES lab website: [Link](https://sees-usf.github.io/)
