#!/usr/bin/env python3
"""Filter a network according to some user specified ."""

import argparse
import sys
import rf_module as rf

def get_args():
    """Get user arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--node-table", dest = "node_table",
                        type = str,
                        help = "Node table as output by direct_network.py")
    parser.add_argument("-e", "--edge-table", dest = "edge_table",
                        type = str,
                        help = "Network as list of edges")
    parser.add_argument("-d", "--d-threshold", dest = "d_threshold",
                        type = float, default = 0.0,
                        help = "Minumum D value of target nodes to include \
                                default = 0.0")
    parser.add_argument("-f", "--f1-score", dest = "f_threshold",
                        type = float, default = 0.9,
                        help = "Minumum F1 score of target nodes to include \
                                default = 0.9")
    parser.add_argument("-o", "--output", dest = "output_file",
                        type = str, help = "Output network file")
    args = parser.parse_args()
    if None in [args.node_table, args.edge_table, args.output_file]:
        parser.print_help(sys.stderr)
        sys.exit(0)
    return [args.node_table, args.edge_table, args.d_threshold,
            args.f_threshold, args.output_file]


def main():
    """
    Filter edges pointing at nodes with D < the threshold or either
    F score < the theshold.
    """
    performance, edge_table, d_threshold, f_threshold, outfile = get_args()
    low_f = {}
    d_table = {}
    for line in rf.get_file_data(performance)[1:]:
        fields = line.split(",")
        low_f[fields[0]] = min(float(fields[16]), float(fields[7]))
        d_table[fields[0]] = float(fields[34])
    edges = rf.get_file_data(edge_table)
    new_edges = [edges[0]]
    for edge in edges[1:]:
        fields = edge.split(",")
        if low_f[fields[0]] > f_threshold and d_table[fields[0]] > d_threshold:
            new_edges.append(edge)
    with open(outfile, "w", encoding = "utf8") as out:
        out.write("\n".join(new_edges) + "\n")

if __name__ == "__main__":
    main()
