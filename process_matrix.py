#!/usr/bin/env python3.6
"""Process the matrix so that it's ready for analysis."""

import sys
import math
import pandas as pd
import rf_module as rf

def get_args():
    """Get user arguments."""
    if len(sys.argv) != 3:
        print("USAGE: python3 process_matrix.py infile outfile\n\n")
        sys.exit()
    return sys.argv[1:]

def write_gene_lists(matrix):
    """
    Write the files - singletons.txt, core_genes.txt, constant_genes.txt.
    """
    constant = list(matrix[matrix.sum(axis = 1) == matrix.shape[1]].index)
    threshold = math.ceil(0.05 * matrix.shape[1])
    core = list(matrix[matrix.sum(axis = 1) >= threshold].index)
    singletons = list(matrix[matrix.sum(axis = 1) == 1].index)
    with open("constant_genes.txt", "w", encoding = "utf-8") as out:
        out.write("\n".join(constant))
    with open("core_genes.txt", "w", encoding = "utf-8") as out:
        out.write("\n".join(core))
    with open("singletons.txt", "w", encoding = "utf-8") as out:
        out.write("\n".join(singletons))

    matrix = matrix.drop(index = constant)
    matrix = matrix.drop(index = singletons)
    return matrix

def collapse_genes(matrix):
    """Collapse genes with the same presence absence pattern."""
    colapsed_rows = pd.DataFrame(columns = matrix.columns)
    identical_sets = {}
    count = 0
    pattern = 0
    indices  = []
    non_unique = matrix[matrix.duplicated(keep = False)]
    non_unique = non_unique.sort_values(by = list(non_unique.columns))
    for i in range(len((matrix[matrix.duplicated(keep = False)].index))):
        if list(non_unique.iloc[i]) == pattern:
            identical_sets["family_group_" + str(count)].append(non_unique.iloc[i].name)
        else:
            colapsed_rows.loc[count] = list(non_unique.iloc[i])
            pattern = list(non_unique.iloc[i])
            count += 1
            indices.append("family_group_" + str(count) + ",,family_group")
            identical_sets["family_group_" + str(count)] = [non_unique.iloc[i].name]

    colapsed_rows.index = indices
    with open("non-unique_genes.csv", "w", encoding = "utf-8") as out:
        for key, value in identical_sets.items():
            out.write(key + "\t" + ",".join(value) + "\n")
    matrix = matrix.drop_duplicates(keep = False)
    matrix = pd.concat([matrix, colapsed_rows])
    return matrix

def collapse_genomes(matrix):
    """Collapse genomes with identical presence absence patterns."""
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
            identical_sets["genome_group_" + str(count)].append(non_unique.iloc[i].name)
        else:
            pattern = list(non_unique.iloc[i])
            colapsed_genomes.loc[count] = list(non_unique.iloc[i])
            count += 1
            indices.append("genome_group_" + str(count))
            identical_sets["genome_group_" + str(count)] = [non_unique.iloc[i].name]

    colapsed_genomes.index = indices
    with open("non-unique_genomes.csv", "w", encoding = "utf-8") as out:
        for key, value in identical_sets.items():
            out.write(key + "\t" + ",".join(value) + "\n")
    matrix = matrix.drop_duplicates(keep = False)
    matrix = pd.concat([matrix, colapsed_genomes])
    matrix = matrix.transpose()
    return matrix


def main():
    """
    Preprocessing Pipeline.

    Remove constant genes - writing their names to a file.
    Remove genes with identical patterns of presence/absence and write to
    a file.
    Remove genomes with identical gene presence/absence patterns and write
    them to a file.
    Write the new matrix to a file.
    """
    infile, outfile = get_args() # pylint: disable=unbalanced-tuple-unpacking
    matrix = pd.read_csv(infile, header = 0, index_col = [0,1,2], dtype = str)
    #remove genes with < 2 present/absent
    matrix = rf.preprocess_df(matrix, 0, 0, 0)
    print("Writing singletons, core and constant genes")
    matrix = write_gene_lists(matrix)

    print("Collapsing identical genes and genomes")
    matrix = collapse_genes(matrix)
    matrix = collapse_genomes(matrix)

    #now write matrix to a file
    print("Writing collapsed matrix")
    matrix.to_csv(outfile, sep = ",", doublequote=False)

if __name__ == "__main__":
    main()
