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

Link to PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/QEST/CRNs/motility_regulation/motility_unb.sm)

Link to JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/QEST/CRNs/motility_regulation/motility_regulation.json) 

The following command runs the framework on this model (from the repo's main directory):

```bash
python  ./Main.py  CRNs/motility_regulation/motility_regulation.json
```