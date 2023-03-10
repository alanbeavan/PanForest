#!/usr/bin/env python3
"""Make an edge table with some extra stats."""

import argparse
import glob
import sys
import pandas as pd
import rf_module as rf

def get_args():
    """Get user arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--matrix", dest = "matrix_file",
                        type = str, help = "Matrix file")
    parser.add_argument("-n", "--network", dest = "network_file",
                        type = str, help = "network file")
    parser.add_argument("-o", "--output", dest = "output_file",
                        type = str, help = "Output table file")
    args = parser.parse_args()
    if None in [args.matrix_file, args.network_file, args.output_file]:
        parser.print_help(sys.stderr)
        sys.exit(0)
    return [args.matrix_file, args.network_file, args.output_file]


def main():
    """
    Write the interactions present in all clusters with extra information.
    each gene in the presence/asbsence of the other.
    end result might be a table w/
    source,target,interaction_type,strength,P(B),P(B|A)P(B|!A),P(A),P(A|B),P(A|!B)P(A+B).
    """
    #for now just do one cluster
    matrix_file, network_file, output_file = get_args()
    matrix = pd.read_csv(matrix_file, header = 0,
                         index_col = [0,1,2], dtype = str)
    matrix = rf.preprocess_df(matrix, 0, 0, 0)
    new_rows = []
    for row in matrix.index:
        new_rows.append(row.split(",")[0])
    matrix.index = new_rows

    interactions = rf.get_file_data(network_file)[1:]

    table_lines = ["Source(A),Target(B),interaction_type,weight,P(B),P(B|A),P(B|!A),P(A),P(A|B),P(A|!B),P(A+B)"]
    matrix = matrix.transpose()
    n_genome = matrix.shape[0]
    for interaction in interactions:
        target, source, interaction_type, weight = interaction.split(",")
        #a_freq = "%.3f"%(len(matrix.loc[matrix[source] == 1])/n_genome)
        a_freq = f'{len(matrix.loc[matrix[source] == 1])/n_genome:.3f}'
        b_freq = f'{len(matrix.loc[matrix[target] == 1])/n_genome:.3f}'
        test_a1 = matrix.loc[matrix[source] == 1]
        test_a0 = matrix.loc[matrix[source] == 0]
#        conditional_b_freq = '%3f'%(sum(test_a1[target])/test_a1.shape[0])
        conditional_b_freq = f'{sum(test_a1[target])/test_a1.shape[0]:.3f}'
        #neg_b_freq = '%3f'%(sum(test_a0[target])/test_a0.shape[0])
        neg_b_freq = f'{sum(test_a0[target])/test_a0.shape[0]:.3f}'
        test_b1 = matrix.loc[matrix[target] == 1]
        test_b0 = matrix.loc[matrix[target] == 0]
#        conditional_a_freq = '%3f'%(sum(test_b1[source])/test_b1.shape[0])
        conditional_a_freq = f'{sum(test_b1[source])/test_b1.shape[0]:.3f}'
        #neg_a_freq = '%3f'%(sum(test_b0[source])/test_b0.shape[0])
        neg_a_freq = f'{sum(test_b0[source])/test_b0.shape[0]:.3f}'
        #p_both = '%3f'%(sum(test_b1[source])/matrix.shape[0])
        p_both = f'{sum(test_b1[source])/matrix.shape[0]:.3f}'
        print(",".join([source,target,interaction_type,weight,
            b_freq,conditional_b_freq,neg_b_freq,
            a_freq,conditional_a_freq,neg_a_freq,p_both]))
        table_lines.append(",".join([source,target,interaction_type,weight,
                                     b_freq,conditional_b_freq,neg_b_freq,
                                     a_freq,conditional_a_freq,neg_a_freq,p_both]))
    with open(output_file, "w", encoding = "utf8") as out:
        out.write("\n".join(table_lines))


if __name__ == "__main__":
    main()
