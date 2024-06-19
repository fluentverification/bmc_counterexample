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


Link to the PRISM model: [Link](https://github.com/fluentverification/bmc_counterexample/blob/QEST/CRNs/enzymatic_futile_cycle/enzym_unb.sm)

Link to the JSON file: [Link](https://github.com/fluentverification/bmc_counterexample/blob/QEST/CRNs/enzymatic_futile_cycle/enzymatic_futile_cycle.json)