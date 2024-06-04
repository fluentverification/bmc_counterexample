### Genetic Circuit0x8E  

Genetic Circuit0x8E is a chemical reaction network of a genetic circuit model consisting of 18 chemical species reacting through 15 reaction channels. This model is defined by the following set of abstracted reactions:

```math
\begin{array}{lll}
    R_1 : \ \textrm{10S}_\textrm{1} \xrightarrow{k_{1}} \emptyset, &
    R_2 : \ \textrm{10S}_\textrm{3} \xrightarrow{k_{1}} \emptyset, \\
    R_3 : \ \textrm{10S}_\textrm{4} \xrightarrow{k_{1}} \emptyset, &
    R_4 : \ \textrm{10S}_\textrm{6} \xrightarrow{k_{1}} \emptyset, \\
    R_5 : \ \textrm{10S}_\textrm{8} \xrightarrow{k_{1}} \emptyset, &
    R_6 : \ \textrm{S}_\textrm{9} + \textrm{S}_\textrm{2} \xrightarrow{k_{2}} \textrm{S}_\textrm{9} + \textrm{S}_\textrm{2} + \textrm{10S}_\textrm{1}, \\
    R_7 : \ \textrm{S}_\textrm{10} + \textrm{S}_\textrm{4} \xrightarrow{k_{3}} \textrm{S}_\textrm{10} + \textrm{S}_\textrm{4} + \textrm{10S}_\textrm{1}, &
    
\end{array}
```
where the initial populations of the species ($S_1$, $S_2$, $S_3$, $S_4$, $S_5$, $S_6$, $S_7$, $S_8$, $S_9$, $S_{10}$, $S_{11}$, $S_{12}$, $S_{13}$, $S_{14}$, $S_{15}$, $S_{16}$, $S_{17}$, $S_{18}$)
are 

```math
\begin{array}{lll}
x_0 = [70, 60, 70, 0, 0, 70, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2].
\end{array}
```
This model uses Hill functions for the propensity of reactions. The reaction rates can be found in the PRISM file.

Link to PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/IEEE/CRNs/circuit0x8E/Circuit0x8E_100to111_unb.sm)

Link to JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/IEEE/CRNs/circuit0x8E/circuit0x8E.json)

The following command runs the framework on this model (from the repo's main directory):

```bash
python  ./Main.py  ./CRNs/circuit0x8E/circuit0x8E.json
```


