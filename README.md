# ReconBench: A Phylogenetic MCMC Benchmark

This is a comprehensive benchmarking software to evaluate phylogenetic MCMC reconstructors such as Historian, BAli-Phy, and others. (Currently only supported for MacOS systems.)

Note: (Rename to recon-bench)

## Overview
Reconstructing ancestral sequence histories is a cornerstone of evolutionary biology, providing insight into the mechanisms of molecular evolution, protein function, and phylogenetic relationships. Accurate ancestral inference requires not only modeling substitutions but also realistically accounting for insertions and deletions (indels), which complicate both alignment and tree estimation. While tools such as Historian have been developed to jointly reconstruct alignments and ancestral states under probabilistic models, systematic benchmarking across diverse evolutionary scenarios remains limited.  
The purpose of this work was to design and implement a comprehensive evaluation framework for Historian, enabling rigorous testing across a wide range of phylogenetic conditions. Using INDELible to simulate evolutionary histories, we generated biologically-derived datasets spanning the four major SCOP protein types containing a diversity of balanced and unbalanced trees, star-like and pectinate structures, variable branch lengths, substitution rates, and indel regimes. These simulations provided a controlled ground truth against which to assess reconstruction accuracy, runtime performance, and robustness.  


## Documentation
If you want to read more about the documentation for ReconBench, a PDF is attached to the repository [].
You can also reach out to Krishna Bhatt [krishbhatt2019@gmail.com] or the Holmes Lab [] for any questions.

## Download

## Requirements
ReconBench requires Python 3.9+. All python dependencies can be installed through pip by

```
pip install -r requirements.txt
```
Additionally, all external tools required for MacOS ARM machines are included within the zipped tools file, which can be extracted by 

```
unzip tools.zip
```

For reference, here is a list of all required software within the 'tools' folder:
* [`ModelTest-NG v0.1.7`](https://github.com/ddarriba/modeltest/releases/tag/v0.1.7)
* [`FastSP v1.7.1`](https://github.com/smirarab/FastSP)
* [`Historian-Mod v0.0.1`](https://github.com/blizzard-labs/historian-mod)
* [`INDELible v1.03`](https://github.com/evolbioinfo/indelible)

Additionally, the following software has to be installed system-wide.

## Running the Evaluation


