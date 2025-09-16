# 2025 Historian & BAli-Phy Evaluation (ReconBench)
This branch contains the results of the 2025 Summer Historian/BAli-Phy. You can find the original ReconBench software [`here`](https://github.com/blizzard-labs/phylo-mcmc-evaluation/tree/2025-summer-evaluation).

## Directory Structure

**`data/matrices`**: Contains all matrices that can be used for the evaluation. The experiments here used the LG08 matrix provided in `data/matrices/lg.json`.

**`data/model_gen`**: Contains the pickled model objects used to generate synthetic datasets for each SCOP type, along with summary statistics and plots available in `data/model_gen/SCOPtypeX/model_info`.

**`data/simulation`**: Contains simulation files and relevant output files from the MCMC software for each of the three experiments. Due to file size constraints, Historian and BAli-Phy output files are compressed into `historian.zip` and `baliphy-1.zip`, respectively. If file sizes are too large, trace files may be seperated into a archive such as `baliphy-trace.zip`.

**`data/results`**: Contains final results tables, including a CSV file including simulated evolution parameters and performance metrics for Historian and BAli-Phy. For experiment 2, there are additional archvies, `historian-plots.zip` and `baliphy-plots.zip` that contains correlation analysis between various, diverse parameters and performance metrics (plots, charts, and more)


## Experimental Summary

### Experiment 1: Representative Sets (Equal Iterations)
For the first experiment, I used four representative parameter sets, one for each respective SCOP type, with indel rates and extension rates equal to those derived in Annabel’s modelfitting of the GGI model. The purpose of this experiment was to analyze Historian and BAli-Phy’s performance on realistic and representative datasets, specifically testing long-run convergence, accuracy, and efficiency over a set number of iterations, independent of any extreme parameters.
* Number of Iterations: 2000 iterations per sequence to fully evaluate convergence
* Number of Sequences: 4 (1 per SCOP dataset) → 8 runs
### Experiment 2: Diverse Sets (Equal Time)
In the second experiment, I used 5 diverse parameters sets per SCOP type, allowing indel rates along with other parameters to vary as per sampling restrictions (discussed previously). Here, the purpose of this experiment was to evaluate how Historian and BAli-Phy’s performs under more varied conditions (indel length, number of indels, sequence length, etc.) Additionally, this limits each program to an equal runtime.
* Maximum Run-time: 2 hours of wall-clock time 
* Number of Sequences: 20 (5 per SCOP dataset) → 40 runs
### Experiment 3: Representative Sets (Equal Time)
Similar to the first experiment, the purpose of this experiment was to analyze Historian and BAli-Phy’s performance on realistic and representative datasets, however over a more realistic constraint of time. This experiment used the same set up as Experiment #1, including four representative parameter sets, one for each respective SCOP type. All parameters are kept the same as well.
* Maximum Run-time: 7 hours of wall-clock time
* Number of Sequences: 4 (1 per SCOP dataset) → 8 runs 

You can find the full details of this evaluation in [`this document`](https://docs.google.com/document/d/1fC3UFOoWkuVDuJioU_jt4Lg4zIOlZc2OwekApCNZfwo/preview).

## Overview of ReconBench
Reconstructing ancestral sequence histories is a cornerstone of evolutionary biology, providing insight into the mechanisms of molecular evolution, protein function, and phylogenetic relationships. Accurate ancestral inference requires not only modeling substitutions but also realistically accounting for insertions and deletions (indels), which complicate both alignment and tree estimation. While tools such as Historian have been developed to jointly reconstruct alignments and ancestral states under probabilistic models, systematic benchmarking across diverse evolutionary scenarios remains limited.  

The purpose of this work was to design and implement a comprehensive evaluation framework for Historian, enabling rigorous testing across a wide range of phylogenetic conditions. Using INDELible to simulate evolutionary histories, we generated biologically-derived datasets spanning the four major SCOP protein types containing a diversity of balanced and unbalanced trees, star-like and pectinate structures, variable branch lengths, substitution rates, and indel regimes. These simulations provided a controlled ground truth against which to assess reconstruction accuracy, runtime performance, and robustness.  

