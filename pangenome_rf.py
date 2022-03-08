#!/usr/bin/env python3.6
"""Train a random forest to identify correlations between genes etc."""

import sys
import os
import math
import random
import argparse
import pandas as pd
import rf_module as rf


def get_args():
    """Get settings from user arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--n-trees", type = int,
                        help = "number of trees in the forest",
                        default = 100, dest = "ntrees")
    parser.add_argument("-d", "--depth", type = int,
                        help = "max depth of trees in the forest",
                        default = 2, dest = "depth")
    parser.add_argument("-m", "--matrix",
                        help = "matrix file", dest = "filename")
    parser.add_argument("-pres", "--min-present", type = float,
                        help = "minimum percentage of genomes featuring a \
                                gene for it to be analysed (5 = 5 percent, \
                                not 0.05)",
                        default = 1.0, dest = "min_present")
    parser.add_argument("-abs", "--min-absent", type = float,
                        help = "minimum percentage of genomes missing a \
                                gene for it to be analysed (5 = 5 percent, \
                                not 0.05)",
                        default = 1.0, dest = "min_absent")
    parser.add_argument("-o", "--output-directory", dest = "output",
                        help = "destiny directory for results")
    parser.add_argument("-r", "--randomise",
                        help = "randomise the gene presence and absence (null\
                                 hypothesis)",
                        action = "store_true")
    parser.add_argument("-t", "--n-threads", dest = "nthreads", default = 1,
                        type = int, help = "number of threads (default - 1)")
    parser.add_argument("-c", "--checkpoint", dest = "checkpoint",
                        default = 0, type = int,
                        help = "continue from checkpoint? provide the number\
                                of genes that have been completed")
    args = parser.parse_args()
    if None in [args.filename, args.output]:
        parser.print_help(sys.stderr)
        sys.exit(0)
    return [args.ntrees, args.depth, args.filename, args.min_present,
            args.min_absent, args.output, args.randomise, args.nthreads,
            args.checkpoint]

def main():
    """
    Pipeline:

    Preprocess matrix:
        Convert the first 3 columns to the row header.
            Currently this returns a warning bc it needs to be a string, not
            a tuple so I could change this
        Convert non-empty entries to 1s and empty to 0.
        Filter the matrix so that only genes that are present in at least
            2 genomes and absent in at least 2 genomes.
    Set up random forest parameters (easy to see at the top of the program).
    Initialise the performance table and the importnace matrices.
    Randomise the order of the genomes (not genes)
        note. I could randomise the genes too to make parallelisation easier
        in future.
    For each gene, perform random forest.
    Split the data into genomes where the gene is present and absent.
        Split each into 75% training and 25% test.
        Separate the y variable (gene presence or absence) for each.
        Fit Classifier using clever other people's code.
        Update performance and importance matrices.
    Add the diagonal to the importance matrices and write to file.
    Update the performance table with various measures.
    """
    ntrees, depth, filename, min_present, min_missing, output,\
            null_h, nthreads, checkpoint = get_args()
    if not os.path.exists(output):
        os.mkdir(output)
    table = pd.read_csv(filename, header = 0, index_col = [0,1,2],
                        dtype = str)
    total_genomes = table.shape[1]
    min_missing = math.ceil(min_missing * total_genomes/100)
    min_present = math.ceil(min_present * total_genomes/100)

    table = rf.preprocess_df(table, null_h, min_missing, min_present)
    imp, performance = rf.init_tables(table)

    #randomise genome order
    n_s = table.shape[1] #number of strains
    table = table[random.sample(list(table.columns), n_s)]
    results = rf.fit_classifiers(table, [imp, performance],
                                 [ntrees, depth, nthreads],
                                 output, checkpoint)
    results[0].to_csv(output + "/imp.csv")
    results[1].to_csv(output + "/performance.csv")

if __name__ == "__main__":
    main()
