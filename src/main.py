import argparse
import os
import sys
import subprocess
import logging
from pathlib import Path
from Bio import Phylo, AlignIO, SeqIO
import utils.general as utils
import model_gen_aa.clean_table
import pandas as pd
from io import StringIO
import time
import string
import random
import re
import concurrent.futures

from main_modelgen import modelConstructor
from main_simulation import evolSimulator
from main_evaluate import main_evaluate


if __name__ == "__main__":
    header = """
=====================================================================================================                                                                                                     
_____ _____ _____ _____    _____                     _               _              _____         _ 
|     |     |     |     |  | __  |___ ___ ___ ___ ___| |_ ___ _ _ ___| |_ ___ ___   |   __|_ _ ___| |
| | | |   --| | | |   --|  |    -| -_|  _| . |   |_ -|  _|  _| | |  _|  _| -_|  _|  |   __| | | .'| |
|_|_|_|_____|_|_|_|_____|  |__|__|___|___|___|_|_|___|_| |_| |___|___|_| |___|_|    |_____|\_/|__,|_|

Phylogenetic MCMC Evaluation Tool v0.1.0, released on 09.08.2025 by the Holmes Lab
Written by Krishna Bhatt (krishbhatt2019@gmail.com)
Latest version: https://github.com/blizzard-labs/phylo-mcmc-evaluation

=====================================================================================================                                                                                                    
    """
    print(header)
    parser = argparse.ArgumentParser(description="Phylogenetic MCMC Evaluation Tool")
    characters = string.ascii_letters + string.digits
    
    help_msg = """
        Phylogenetic MCMC Evaluation Tool
        Usage: python src/main.py --mode <mode> --input <input_file> --output <output_file> [OPTIONS]
        
        
        
        Step-by-step guide can be found on GitHub: https://github.com/blizzard-labs/phylo-mcmc-evaluation 
    """
    
    parser.add_argument("--mode", type=str, choices=["modelgen", "simulate", "evaluate"],
                        help="Mode of operation: 'modelgen' to generate models, 'simulate' to run simulations, 'evaluate' to evaluate results",
                        required=True)
    parser.add_argument("--actions", type=str, nargs='+',
                        help="actions to perform in the selected mode",
                        required=False, default=[])
    parser.add_argument("--input", type=str,
                        help="Input file or folder path",
                        required=False)
    parser.add_argument("--label", type=str,
                        help="Label for the current run, used in naming output files",
                        required=False, default=''.join(random.choice(characters) for _ in range(7)))
    parser.add_argument("--help", action="help", 
                        help="Show this help message and exit")
    
    
    args = parser.parse_args()
    print(f"Selected mode: <{args.mode}>")
    
    #*===========================================================================================
    #? Modelgen Mode Functionality
    #*===========================================================================================
    
    if args.mode == "modelgen":        
        supported_actions = ['cleanup-pfam', 'extract-subst-params', 'extract-top-params',
                             'cleanup-params', 'generate-model', 'sample-model']
        actions = [action.lower() in supported_actions for action in args.actions]
        
        if not all(actions):
            print(f"Error: Unsupported actions found in {args.actions}. Supported actions are: {supported_actions}")
            sys.exit(1)
        
        if len(args.actions) == 0:
            actions = supported_actions.copy()
            
        start = time.time()
        #! Modeltest-NG can be run on linux, however to maintainc consistency, this program is based on MacOS
        
        try:
            mc = modelConstructor('osx', args.label, args.input, 
                                args.input.replace('alignments', 'protein_evolution_parameters.csv'), log=False)
        except Exception as e:
            print(f"Error initializing modelConstructor: {e}\nDouble-check the input paths/format")
            sys.exit(1)
        
        for action in args.actions:
            print(f"Running action: <{action}>...")
            if action == 'cleanup-pfam':
                try:
                    mc.cleanup_pfam_data()
                    print(f'Cleaned data- ELAPSED TIME: {time.time() - start}============================')
                except Exception as e:
                    print(f"Error during cleanup_pfam_data: {e}\nDouble-check the input paths/format")
                    sys.exit(1)
            
            elif action == 'extract-subst-params':
                try:
                    mc.extract_substitution_params()
                    print(f'Extracted substitution params- ELAPSED TIME: {time.time() - start}============================')
                except Exception as e:
                    print(f"Error during extract_substitution_params: {e}\nDouble-check the input paths/format")
                    sys.exit(1)
            
            elif action == 'extract-top-params':
                try:
                    mc.extract_top_params()
                    print(f'Extracted topology params- ELAPSED TIME: {time.time() - start}============================')
                except Exception as e:
                    print(f"Error during extract_top_params: {e}\nDouble-check the input paths/format")
                    sys.exit(1)
            
            elif action == 'cleanup-params':
                try:
                    mc.cleanup_params(mc.params_file)
                    print(f'Cleaned parameters file- ELAPSED TIME: {time.time() - start}============================')
                except Exception as e:
                    print(f"Error during cleanup_params: {e}\nDouble-check the input paths/format")
                    sys.exit(1)
            
            elif action == 'generate-model':
                try:
                    mc.generate_model(mc.params_wc_file)
                    print(f'Generated model- ELAPSED TIME: {time.time() - start}============================')
                except Exception as e:
                    print(f"Error during generate_model: {e}\nDouble-check the input paths/format")
                    sys.exit(1)
            
            elif action == 'sample-model':
                try:
                    mc.sample_model(n_samples=1000)
                    print(f'Sampled model- ELAPSED TIME: {time.time() - start}============================')
                except Exception as e:
                    print(f"Error during sample_model: {e}\nDouble-check the input paths/format")
                    sys.exit(1)
            
            print(f"Completed action: <{action}>.")
        print(f"All actions completed. Total ELAPSED TIME: {time.time() - start}============================")
    
    #*===========================================================================================
    #? Simulate Mode Functionality
    #*===========================================================================================
    
    elif args.mode == "simulate":
        supported_actions = ['generate_tree_topologies', 'simulate_evolution', 'cleanup_folders',
                             'run_mcmc']
        
        actions = [action.lower() in supported_actions for action in args.actions]
        
        if not all(actions):
            print(f"Error: Unsupported actions found in {args.actions}. Supported actions are: {supported_actions}")
            sys.exit(1)
        
        if len(args.actions) == 0:
            actions = supported_actions.copy()
            
        start = time.time()
        
        try:
            es = evolSimulator(args.input, tag=args.label, consensus_tree_file='none')
        except Exception as e:
            print(f"Error initializing evolSimulator: {e}\nDouble-check the input paths/format")
            sys.exit(1)
        
        for action in args.actions:
            print(f"Running action: <{action}>...")
            if action == 'generate_tree_topologies':
                try:
                    es.generate_treetop_with_params(max_iterations=1000)
                    print(f'Generated tree topologies- ELAPSED TIME: {time.time() - start}============================')
                except Exception as e:
                    print(f"Error during generate_tree_topologies: {e}\nDouble-check the input paths/format")
                    sys.exit(1)
            
            elif action == 'simulate_evolution':
                try:
                    es.runIndelible()
                    print(f'Simulated evolution- ELAPSED TIME: {time.time() - start}============================')
                except Exception as e:
                    print(f"Error during simulate_evolution: {e}\nDouble-check the input paths/format")
                    sys.exit(1)
            
            elif action == 'cleanup_folders':
                try:
                    es.cleanupSimFolders()
                    print(f'Cleaned up folders- ELAPSED TIME: {time.time() - start}============================')
                except Exception as e:
                    print(f"Error during cleanup_folders: {e}\nDouble-check the input paths/format")
                    sys.exit(1)
            
            elif action == 'run_mcmc':
                try:
                    es.runBenchmark()
                    print(f'Ran MCMC analyses- ELAPSED TIME: {time.time() - start}============================')
                except Exception as e:
                    print(f"Error during run_mcmc: {e}\nDouble-check the input paths/format")
                    sys.exit(1)
            
            print(f"Completed action: <{action}>.")
        print(f"All actions completed. Total ELAPSED TIME: {time.time() - start}============================")
        
    #*===========================================================================================
    #? Evaluate Mode Functionality
    #*===========================================================================================
    
    elif args.mode == "evaluate":
        supported_actions = ['trace_parser', 'clean_treestat', 'convergence',
                               'comparison', 'summarize']
        
        actions = [action.lower() in supported_actions for action in args.actions]
        
        if not all(actions):
            print(f"Error: Unsupported actions found in {args.actions}. Supported actions are: {supported_actions}")
            sys.exit(1)
        
        if len(args.actions) == 0:
            actions = supported_actions.copy()
            
        start = time.time()
        
        for action in args.actions:
            print(f"Running action: <{action}>...")
            try:
                main_evaluate(args.input, action)
                print(f'Completed action: <{action}>- ELAPSED TIME: {time.time() - start}============================')
            except Exception as e:
                print(f"Error during {action}: {e}\nDouble-check the input paths/format")
                sys.exit(1)
        print(f"All actions completed. Total ELAPSED TIME: {time.time() - start}============================")