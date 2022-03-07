#!/usr/bin/env python3.6
"""Calculate if a relationship is negative or positive for all edges."""

import sys
import pandas as pd
import rf_module as rf

def get_args():
    """Get user arguments."""
    if len(sys.argv) != 4:
        print("USAGE: python3 direct_network.py in_graph, pres_abs_matrix,\
                outfile")
        sys.exit()
    return sys.argv[1:]

def get_file_data(filename):
    """Stores the lines of the program with name filename as a list."""
    with open(filename, encoding = "utf8") as in_file:
        lines = []
        for line in in_file:
            lines.append(line.rstrip("\n"))
    return lines

def calculate_proportions(matrix):
    """Calculate the background abundance of all genes."""
    props = {}
    n_strain = matrix.shape[1]
    counts = matrix.sum(axis = 1)
    for index in counts.index:
        props[index] = counts.loc[index]/n_strain
    print(props)
    return props

def main():
    """
    Write a new network file with pp for +ve intecations and nn for -ve.

    Read in the presence absence matrix and the network file.
    Calculate the proportion of genomes with each gene regardless of other
        genes, saving as a dictionary.
    For each line of the network file, calculate the proportion of genomes
        with the target gene present when the source gene is also present.
    Compare this with the background rate. If it is higher, the interavtion
        is positive so the line of the network can be added to the new lines
        as is. Otherwise, it is a negative interaction so "pp" should be
        replaced with nn.
    Write the new network file.
    """
    # pylint: disable=unbalanced-tuple-unpacking
    # pylint: disable=consider-using-enumerate
    infile, matrix_file, outfile = get_args()
    matrix = pd.read_csv(matrix_file, header = 0,
                         index_col = [0,1,2], dtype = str)
    matrix = rf.preprocess_df(matrix, 0, 0, 0)
    index_dict = {}
    for i in range(len(matrix.index)):
        index_dict[matrix.index[i]] = matrix.index[i].split(",")[0]
    matrix = matrix.rename(index = index_dict)
    props = calculate_proportions(matrix)
    matrix = matrix.transpose()

    network = get_file_data(infile)
    newlines = []
    for line in network[1:]:
        target, source = line.split(",")[0:2]
        weight = line.split(",")[3]
        test = matrix.loc[matrix[source] == 1]
        if sum(test[target])/test.shape[0] < props[target]:
            newlines.append(",".join([target,source,"nn",weight]))
        else:
            newlines.append(line)
    with open(outfile, "w", encoding = "utf8") as out:
        out.write("Target,Source,InteractionType,Weight\n")
        out.write("\n".join(newlines))


if __name__ == "__main__":
    main()
