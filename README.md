# ReconBench: A Phylogenetic MCMC Benchmark

This is a comprehensive benchmarking software to evaluate phylogenetic MCMC reconstructors such as Historian, BAli-Phy, and others. (Currently only supported for MacOS systems.)

## Overview
Reconstructing ancestral sequence histories is a cornerstone of evolutionary biology, providing insight into the mechanisms of molecular evolution, protein function, and phylogenetic relationships. Accurate ancestral inference requires not only modeling substitutions but also realistically accounting for insertions and deletions (indels), which complicate both alignment and tree estimation. While tools such as Historian have been developed to jointly reconstruct alignments and ancestral states under probabilistic models, systematic benchmarking across diverse evolutionary scenarios remains limited.  

The purpose of this work was to design and implement a comprehensive evaluation framework for Historian, enabling rigorous testing across a wide range of phylogenetic conditions. Using INDELible to simulate evolutionary histories, we generated biologically-derived datasets spanning the four major SCOP protein types containing a diversity of balanced and unbalanced trees, star-like and pectinate structures, variable branch lengths, substitution rates, and indel regimes. These simulations provided a controlled ground truth against which to assess reconstruction accuracy, runtime performance, and robustness.  


## Documentation
If you want to read more about the documentation for ReconBench, a PDF is attached to [`the repository`](https://github.com/blizzard-labs/phylo-mcmc-evaluation/blob/main/docs/2025_report.pdf), or viewable in [`Google Drive`](https://docs.google.com/document/d/1fC3UFOoWkuVDuJioU_jt4Lg4zIOlZc2OwekApCNZfwo/preview).
You can also reach out to Krishna Bhatt [krishbhatt2019@gmail.com] or the Holmes Lab [holmeslab.org] for any questions.
* [`2025 Summer Evaluation Results`](https://github.com/blizzard-labs/phylo-mcmc-evaluation/tree/2025-summer-evaluation): Comparison of Historian and BAli-Phy conducted during the summer of 2025. All details on this initial evaluation can be found in the attached PDF above.

## Download
You can download Reconbench with 
```
git clone https://github.com/blizzard-labs/phylo-mcmc-evaluation.git
```
After installing all requirements, you can run the script with `python src/main.py --[OPTIONS]`.

## Requirements
ReconBench requires Python 3.9+. All python dependencies can be installed through pip by

```
pip install -r requirements.txt
```
Additionally, all external tools required for MacOS ARM machines are included within the zipped tools file, which can be extracted by 

```
unzip tools.zip
```

For reference, here is a list of all required software within the 'tools' folder (ensure naming matches format in tools.zip):
* [`ModelTest-NG v0.1.7`](https://github.com/ddarriba/modeltest/releases/tag/v0.1.7)
* [`FastSP v1.7.1`](https://github.com/smirarab/FastSP)
* [`Historian-Mod v0.0.1`](https://github.com/blizzard-labs/historian-mod)
* [`INDELible v1.03`](https://github.com/evolbioinfo/indelible)

Additionally, the following software has to be installed and accessible via your terminal `$PATH` variable:
* [`BAli-Phy v4.1`](https://www.bali-phy.org/download.php)
* [`BEAST v2.7.8`](https://www.beast2.org/)

After installing the BEAST software, open BEAST's "AppLauncher" app and install TreeStat2.

## Running the Evaluation

ReconBench has three major functionalities, that can be accessed by three modes:
* `modelgen`: Used to generate an bio-informed model for information. Given alignment and tree data for a phylogenetic group, generates a multivariate distribution of various evolutionary parameters (sequence, topology, etc.) This can be sampled to generate highly realistic synthetic datasets.
* `simulate`: Used to generate sequences and run MCMC. Given datasets of parameters, generates sequences (with INDELible) and trees (through simulated annealing) to fit the data. This data is organized and fed to the MCMC software, recording trace information and wall-clock time.
* `evaluate`: Used to generate final summary statistics for each software. Given output files from Historian and BAli-Phy, parses traces and computes MCMC statistics for topology and sequence. Additionally computes the posterior decoding alignment and CCD-1 MAP tree to compare RF, SP, and more scores to the ground truth. Organizes all data and analyzes correlations.

Here is the basic structure of a call to ReconBench:
```
python src/main.py --mode <mode> --actions <action1 action2 ...> --input <filepath> --label <output_header> [--OPTIONS]
```

## Step-by-step Guide

This short tutorial will give you an example of how to use ReconBench using the example dataset found in `data/model_gen/example`. This was a set of 10 alignments collected from Pfam and their respective trees (generated with Pfam's FastTree program) organized into two folders: `alignments/` and `trees/`. Any input to ReconBench must be provided in the same organization. 

Specifically for Pfam data, we provide a cleaning function to convert to standard FASTA and Newick format. You can call this through
```
python src/main.py --mode modelgen --input data/model_gen/example/alignments --label example --actions cleanup-pfam
```
Note that for the `modelgen` mode, the alignment folder path is specified under the `--input` flag, the dataset's name is given by the `--label` flag, and `--actions` specifies the actions to be taken.

We can then continue to extract parameters to generate and sample from a model, all in one command.

```
python src/main.py --mode modelgen --input data/model_gen/example/alignments --label example --actions extract-subst-params extract-top-params cleanup-params generate-model sample-model
```

Note that multiple actions can be placed together, seperated by spaces. Additionally note that the model cannot actually be generated due to the low number of alignments/trees in our example dataset. For now, I have selected parameters in the `data/model_gen/example/experiment_parameters.csv`. Ensure that after you have sampled your parameters, the match the same format as this file for the next step of simulation. If you have a pickled models object that you want to load in, ensure it is named as `data/model_gen/example/bd_models.pkl`.

**WARNING: Extracting substitution or topology parameters may fail for certain trees/alignments.** This is a known issue- for now, remove all failed sets from the original dataset and try again.

Now, we can move onto simulating datasets. Similarly, you can generate trees, simulate evolution, and clean up your data for all of your experimental parameter sets in one command:

```
python src/main.py --mode simulate --input data/model_gen/example/experiment_parameters.csv --label example --actions generate-tree-topologies simulate-evolution cleanup-folders
```

Note here that the `--input` flag for the simulate mode specifies the path to the CSV file containing all of the parameter datasets. We can now run Historian/BAli-Phy on the ground truths, specifying the maximum run time and maximum number of iterations (per sequence) with the `--timeout` and `--max_iter` flags respectively. The programs will die upon reaching whichever criteria occurs first.

```
python src/main.py --mode simulate --input data/model_gen/example/experiment_parameters.csv --label example --actions run-mcmc --timeout 0.083  --max_iter 30
```

After running the software, we can parse the output trace files with

```
python src/main.py --mode evaluate --input data/simulation --actions trace_parser
```

Note that no `--label` flag is required here and the `--input` flag specifies the location of the simulation folder in this mode. After parsing the traces, open the TreeStat2 GUI software through the BEAST AppLauncher. Here, you can select the following parameters: `CCD1 RF distance`, `CCD1 information content (log(p))`, `Colless tree-imbalance`, `Tree Length`, `Tree Topology`. Then, for each sequence, upload the Historian and BAli-Phy's tree trace named `parsed_trace.log.trees` and `cleaned.trees`, located in their respective folders for each sequence in the simulation folder. Save the TreeStat2 output to the same folder (`historian`/`baliphy-1`), titled `treetrace.log`. This should produce a CCD-1 MAP tree as well.

After running TreeStat on all of the sequences, run the following to complete the evaluation:
```
python src/main.py --mode evaluate --input data/simulation --actions clean_treestat convergence comparison summarize
```

You can find the output of an fully run example with real datasets [`here`](https://github.com/blizzard-labs/phylo-mcmc-evaluation/tree/2025-summer-evaluation).