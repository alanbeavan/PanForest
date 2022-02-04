#!/usr/bin/env python3.6
"""Train a random forest to identify correlations between genes etc."""

from rf_module import *
import pandas as pd
import numpy as np
import sklearn
import random
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
import my_module as mod

def settings():
    """Define the settings. Modify as you want."""
    ntrees = 100
    depth = 3
    filename = "gene_presence_absence.csv"
    null_h = False
    min_missing = 10
    min_present = 10
    return ntrees, depth, filename, null_h, min_missing, min_present

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
    Maybe plot and save some histograms and/or boxplots.                        X
    """
    ntrees, depth, filename, null_h, min_missing, min_present = settings()
    df = pd.read_csv(filename, header = 0, index_col = [0,1,2], dtype = str)
    df = preprocessDf(df, null_h, min_missing, min_present)
    ns = df.shape[1] #number of strains
    ng = df.shape[0] #number of genes (gene families)
    imp = pd.DataFrame(0.0, index = np.arange(ng), columns = df.index.values)
    imp.index = df.index.values
    metrics = ['count', 'TPte', 'FPte', 'FNte', 'TNte', 'Ete', 'Ate', 'P1te',
               'P0te', 'Pte', 'R1te', 'R0te', 'Rte', 'F1te', 'F0te', 'Fte',
               'TPtr','FPtr', 'FNtr', 'TNtr', 'Etr', 'Atr', 'P1tr', 'P0tr',
               'Ptr', 'R1tr', 'R0tr', 'Rtr', 'F1tr', 'F0tr', 'Ftr']
    performance = pd.DataFrame(0.0, index = np.arange(ng), columns = metrics)
    performance.index = df.index.values
    df = df[random.sample(list(df.columns), ns)] #randomise genome order
    fit_classifiers(df, imp, performance, ns, ng, ntrees, depth)
    imp.to_csv("imp.csv")
    performance.to_csv("performance.csv")

if __name__ == "__main__":
    main()
