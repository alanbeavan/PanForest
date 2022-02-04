#!/usr/bin/env python3.6
"""Process the matrix so that it's ready for analysis."""

from rf_module import *
import math
import pandas as pd
import my_module as mod

def main():
    """
    Preprocessing Pipeline.
    
    Remove constant genes - writing their names to a file.
    Remove genes with identical patterns of presence/absence and write to
    a file.
    Remove genomes with identical gene presence/absence patterns and write
    them, you guessed it, to a file.
    Write the new matrix to a novel.

    Of course I mean a file.
    """
    matrix = pd.read_csv("all_merged/gene_presence_absence.csv", header = 0,\
            index_col = [0,1,2], dtype = str)
    #remove genes with < 2 present/absent
    matrix = preprocessDf(matrix, 0, 0, 0)
    constant = list(matrix[matrix.sum(axis = 1) == matrix.shape[1]].index)
    threshold = math.ceil(0.05 * matrix.shape[1])
    core = list(matrix[matrix.sum(axis = 1) >= threshold].index)
    singletons = list(matrix[matrix.sum(axis = 1) == 1].index)
    out = open("constant_genes.txt", "w")
    out.write("\n".join(constant))
    out.close()
    out = open("core_genes.txt", "w")
    out.write("\n".join(core))
    out.close()
    out = open("singletons.txt", "w")
    out.write("\n".join(singletons))
    out.close()

    matrix = matrix.drop(index = constant)
    matrix = matrix.drop(index = singletons)

    colapsed_rows = pd.DataFrame(columns = matrix.columns)
    identical_sets = {}
    count = 0
    pattern = 0
    indices  = []
    non_unique = matrix[matrix.duplicated(keep = False)]
    
    #print("blurt")
    #print(len(non_unique.index))
    #while len(non_unique.index) >= 1:
    #    print("count")
    #    print(count)
    #    pattern = list(non_unique.iloc[0])
    #    genes = list(non_unique[non_unique == pattern].dropna().index)
    #    identical_sets["group_" + str(count)] = genes
    #    indices.append("group_" + str(count))
    #    count += 1
    #    non_unique = non_unique.drop(genes, axis = 0)
    #    print("remainder")
    #    print(len(non_unique.index))

    non_unique = non_unique.sort_values(by = list(non_unique.columns))
    for i in range(len((matrix[matrix.duplicated(keep = False)].index))):        
        if list(non_unique.iloc[i]) == pattern:
            identical_sets["group_" + str(count)].append(non_unique.iloc[i].name)
        else:
            colapsed_rows.loc[count] = list(non_unique.iloc[i])
            pattern = list(non_unique.iloc[i])
            count += 1
            indices.append("group_" + str(count))
            identical_sets["group_" + str(count)] = [non_unique.iloc[i].name]
            print(colapsed_rows)

    colapsed_rows.index = indices
    out = open("non-unique_genes.csv", "w")
    for key, value in identical_sets.items():
        out.write(key + "\t" + ",".join(value) + "\n")
    out.close()
    matrix = matrix.drop_duplicates(keep = False)
    matrix = pd.concat([matrix, colapsed_rows])

    #now transpose and remove duplicate genomes
    matrix = matrix.transpose()
    colapsed_genomes = pd.DataFrame(columns = matrix.columns)
    identical_sets = {}
    count = 0
    pattern = 0
    indices  = []
    non_unique = matrix[matrix.duplicated(keep = False)]
    non_unique = non_unique.sort_values(by = list(non_unique.columns))
    
    for i in range(len((matrix[matrix.duplicated(keep = False)].index))):
        if list(non_unique.iloc[i]) == pattern:
            identical_sets["g_group_" + str(count)].append(non_unique.iloc[i].name)
        else:
            pattern = list(non_unique.iloc[i])
            colapsed_genomes.loc[count] = list(non_unique.iloc[i])
            count += 1
            indices.append("g_group_" + str(count))
            identical_sets["g_group_" + str(count)] = [non_unique.iloc[i].name]
            print(colapsed_genomes)
    
    colapsed_genomes.index = indices
    out = open("non-unique_genomes.csv", "w")
    for key, value in identical_sets.items():
        out.write(key + "\t" + ",".join(value) + "\n")
    out.close()
    matrix = matrix.drop_duplicates(keep = False)
    matrix = pd.concat([matrix, colapsed_genomes])
    matrix = matrix.transpose()

    #now write matrix to a file
    matrix.to_csv("collapsed_matrix.csv", sep = ",")

if __name__ == "__main__":
    main()
