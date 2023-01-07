# fuzzy-rainbow

This [python tool](https://github.com/secgroup/fuzzy-rainbow/blob/main/calc.py) computes the parameter for the reference design proposed in the paper:

Leonardo Veronese, Francesco Palmarini, Riccardo Focardi and Flaminia L. Luccio. _A fast and cost-effective design for FPGA-based fuzzy rainbow tradeoffs_.
Proceedings of the [8th International Conference on Information Systems Security and Privacy (ICISSP'22)](https://icissp.scitevents.org/?y=2022), 9-11 February, 2022.

The tool can be used for:

- computing the optimal parameter set for a given hardware and attack constraints;
- estimating the precomputation time and online performance;
- reproducing and validating the results presented in the paper.

## Reproducing the results in the paper

To reproduce the parameters and the estimated performance for DES and A5/1 presented in the paper set the `CIPHER` variable to one of the following values
and execute the tool:

```python
CIPHER = 'DES'       # (Table 1)
CIPHER = 'A5/1 90%'  # (Table 2, first column)
CIPHER = 'A5/1 99%'  # (Table 2, second column)
```
## General usage

To use the tool for a specific application:

1. define a new set of global parameters with hardware and attack constraints;
2. tune compression parameters by adjusting the values based on memory constraints and the tool output.

**NOTE:** the tool automatically finds the optimal `s` value starting from value `15` and iterating the search until a fixpoint is reached.
