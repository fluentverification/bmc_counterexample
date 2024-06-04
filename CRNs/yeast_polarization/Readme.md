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