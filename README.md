
# Stochastic Analysis of Infinite-State Chemical Reaction Networks

The framework takes PRISM models of CRNs and a CSL property of the form $P=? [\textrm{true} \; \textrm{U}<=\textrm{T} \;  (\textrm{X}=\theta)]$ as input, and generates and iteratively increasing lower-bound for the probability of the property on the CRN model.
At each iteration, the lower-bound probabilty and the size of the state-space on which this probability is computed is printed to the console.

# Installation

- Prerequisites

	- Stormpy

		- Stormpy is a set of python bindings for the STORM’s C++ API. The framework uses this library for parsing models. Stormpy can be installed following the instructions showed [here](https://moves-rwth.github.io/stormpy/installation.html#installation-steps).

	- PRISM

		- PRISM is used as the back-end probabilistic model-checker. The location to the PRISM binary should be passed as a parameter to the framework. Details on installation of PRISM can be found [here](https://www.prismmodelchecker.org/manual/InstallingPRISM/Instructions).

	- z3-py

		- This framework uses z3 and its python bindings z3-py as the underlying SMT-solver. z3-py can added to the python interpreter by running the following command: `pip install z3-solver` . More details on installation of z3-py can be found [here](https://github.com/Z3Prover/z3?tab=readme-ov-file#z3-bindings).

- **Note: Both Stormpy and z3-py should be installed in the same python virtual environment.**

- Cloning the python scripts
	
	- After prerequisite libraries are installed, the python scripts can be cloned and accessed by following commands:
		 ```bash
		git clone --depth 1 --branch IEEE https://github.com/fluentverification/bmc_counterexample.git
		cd bmc_counterexample
		```


# Running

The framework accepts CTMC models written in [PRISM modeling language](https://www.prismmodelchecker.org/manual/ThePRISMLanguage/Introduction) as input. Currently, the framework only accepts a restricted class of models where:

1. All commands in the model are labelled.

2. There is no module-renaming in the model.

3. Each command in the model must have exactly one update.

To run the framework, a JSON file containing the necessary parameters should be passed to the python interpreter. This JSON file should follow this format:

```json
{
"model_name" : (string)  arbitrary  name  given  to  the  model,
"model_path" : (string)  the  path  to  the  PRISM  model,
"starting_bound" : (string) length of shortest witness trace. Can be set to 1 if not known,
"prism_binary" : (string) path to the PRISM binary,
"csl_property" : (string)  CSL  property  that  is  going  to  be  checked,
"target_variable" : (string)  the  variable  of  interest  in  the  CSL  property,
"target_value" : (string)  the  value  of  the  variable  of  interest  in  the  CSL  property,
"model_check_step" : (string) only necessary for XBF. PRISM is going to be called after the size of state space has increased by at least this value,
"mode" : (string) one of 1)bdfs_1 2)bdfs_1_plus 3)bdfs_all
}
```

To run the framework, this JSON file should be simply passed as an argument to the [Main.py](http://Main.py) script:

```bash
python  Main.py  file.json
```  

# Case Studies  

### Enzymatic Futile Cycle

Enzymatic Futile Cycle is a chemical reaction network consisting of 6 chemical species reacting through 6 reaction channels. The model is defined by the following set of reactions:

```math
\begin{array}{lll}
    R_1 : \ \textrm{S}_1 + \textrm{S}_2 \xrightarrow{1.0} \textrm{S}_3,~~~ &
    R_2 : \ \textrm{S}_3 \xrightarrow{1.0} \textrm{S}_1 + \textrm{S}_2,\\
    R_3 : \ \textrm{S}_3 \xrightarrow{0.1} \textrm{S}_1 + \textrm{S}_5, ~~~ &
    R_4 : \ \textrm{S}_4 + \textrm{S}_5 \xrightarrow{1.0} \textrm{S}_6,\\
    R_5 : \ \textrm{S}_6 \xrightarrow{1.0} \textrm{S}_4 + \textrm{S}_5,~~~ &
    R_6 : \ \textrm{S}_6 \xrightarrow{0.1} \textrm{S}_4 + \textrm{S}_2
\end{array}
```
where the initial populations of species $(S_1, S_2, S_3, S_4, S_5, S_6)$ are 
```math
\begin{array}{lll}
    x_0 = [1, 50, 0, 1, 50, 0].
\end{array}
```


Link to the PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/IEEE/CRNs/enzymatic_futile_cycle/enzym_unb.sm)

Link to the JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/IEEE/CRNs/enzymatic_futile_cycle/enzymatic_futile_cycle.json)

The following command runs the framework on this model (from the repo's main directory):

```bash
python  ./Main.py  ./CRNs/enzymatic_futile_cycle/enzymatic_futile_cycle.json
```

### Motility Regulation

Motility regulation is a chemical reaction network of a gene regulatory network consisting of 9 chemical species reacting through 12 reaction channels. The model is defined by the following set of reactions:

```math
\begin{array}{lll}
    R_1 : \ \textrm{codY} \xrightarrow{0.1} \textrm{codY} + \textrm{CodY}, &
    R_2 : \ \textrm{CodY} \xrightarrow{0.0002} \emptyset, \\
    R_3 : \ \textrm{flache} \xrightarrow{1} \textrm{flache} + \textrm{SigD}, &
    R_4 : \ \textrm{SigD} \xrightarrow{0.0002} \emptyset, \\
    R_5 : \ \textrm{SigD\_hag} \xrightarrow{1} \textrm{SigD} + \textrm{hag} + \textrm{Hag}, &
    R_6 : \ \textrm{Hag} \xrightarrow{0.0002} \emptyset, \\
    R_7 : \ \textrm{SigD} + \textrm{hag} \xrightarrow{0.01} \textrm{SigD\_hag}, &
    R_8 : \ \textrm{SigD\_hag} \xrightarrow{0.1} \textrm{SigD} + \textrm{hag}, \\
    R_9 : \ \textrm{CodY} + \textrm{flache} \xrightarrow{0.02} \textrm{CodY\_flache},&
    R_{10} : \ \textrm{CodY\_flache} \xrightarrow{0.1} \textrm{CodY} + \textrm{flache}, \\
    R_{11} : \ \textrm{CodY} + \textrm{hag} \xrightarrow{0.01} \textrm{CodY\_hag}, &
    R_{12} : \ \textrm{CodY\_hag} \xrightarrow{0.1} \textrm{CodY} + \textrm{hag} 
\end{array}
```

where the initial populations of the species 
(codY, CodY, flache, SigD, SigD\_hag, hag, Hag, CodY\_flache, CodY\_hag)
are 

```math
\begin{array}{lll}
    x_0 = [1, 10, 1, 10, 1, 1, 10, 1, 1].
\end{array}
```

Link to PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/IEEE/CRNs/motility_regulation/motility_unb.sm)

Link to JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/IEEE/CRNs/motility_regulation/motility_regulation.json) 

The following command runs the framework on this model (from the repo's main directory):

```bash
python  ./Main.py  CRNs/motility_regulation/motility_regulation.json
```

### Yeast Polarization

Yeast Polarization is a chemical reaction network model consisting of 7 chemical species reacting through 8 reaction channels. The model is defined by the following set of reactions:

```math
\begin{array}{ll}
    R_1 : \ \emptyset \xrightarrow{0.0038} \textrm{R}, &
    R_2 : \ \textrm{R} \xrightarrow{4.00\times 10^{-4}} \emptyset,\\
    R_3 : \ \textrm{L} + \textrm{R} \xrightarrow{0.042} \textrm{RL} + \textrm{L}, &
    R_4 : \ \textrm{RL} \xrightarrow{0.0100} \textrm{R},\\
    R_5 : \ \textrm{RL} + \textrm{G} \xrightarrow{0.011} \textrm{G}_\textrm{a} + \textrm{G}_{\textrm{bg}},&
    R_6 : \ \textrm{G}_\textrm{a} \xrightarrow{0.100} \textrm{G}_\textrm{d}, \\
    R_7 : \ \textrm{G}_\textrm{d} + \textrm{G}_{\textrm{bg}} \xrightarrow{1.05\times 10^{3}} \textrm{G},&
    R_8 : \ \emptyset \xrightarrow{3.21} \textrm{RL} 
\end{array}
```
where the initial populations of species $(R, L, RL, G, G_{a}, G_{bg}, G_d)$ are 

```math
\begin{array}{lll}
x_0 = [50, 2, 0, 50, 0, 0, 0].
\end{array}
```

Link to PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/IEEE/CRNs/yeast_polarization/yeast_unb.sm) 

Link to JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/IEEE/CRNs/yeast_polarization/yeast_polarization.json)

The following command runs the framework on this model(from the repo's main directory):

```bash
python  Main.py  ./CRNs/yeast_polarization/yeast_polarization.json
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