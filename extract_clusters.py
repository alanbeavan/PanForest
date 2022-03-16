#!/usr/bin/env python3.6
"""Cluster a network and write the clusters to a file."""

import sys
import argparse
import os
import pandas as pd
import networkx as nx
from networkx.algorithms.community import louvain_communities
from networkx.algorithms.community import greedy_modularity_communities
from networkx.algorithms.community import asyn_lpa_communities

def get_args():
    """Get settings from user arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str,
                        help = "Input network file", dest = "infile")
    parser.add_argument("-e", "--edge-type", type = str,
                        help = "Interaction type (pp, nn). \
                                Default = null (both included)",
                        default = False, dest = "edge_type")
    parser.add_argument("-m", "--method", type = str, dest = "method",
                        help = "Method to cluster graph (modularity (greedy),\
                                 louvain (iterative), default = modularity",
                        default = "modularity")
    parser.add_argument("-o", "--output", type = str, default = "clusters",
                        help = "Output directory, default = clusters/",
                        dest = "outdir")
    args = parser.parse_args()
    if None in [args.infile]:
        parser.print_help(sys.stderr)
        sys.exit(0)
    return [args.infile, args.edge_type, args.method, args.outdir]


def main():
    """
    Read in edge table
    Cluster graph
    Write clusters to a file
    """
    infile, edge_type, method, outdir = get_args()
    edges = pd.read_csv(infile, header = [0])
    #edges = edges.head(100)
    if edge_type:
        edges = edges[edges["InteractionType"] == edge_type]
    network = nx.from_pandas_edgelist(edges, source = "Source",
                                      target = "Target",
                                      edge_attr = ["InteractionType", "Weight"],
                                      create_using=nx.DiGraph())
    if method == "louvain":
        clusters = louvain_communities(network, weight = "Weight",
                                       threshold = 0.0001)
    elif method == "modularity":
        clusters = greedy_modularity_communities(network, weight = "Weight")
    else:
        print("method must be louvain or modularity. Try again")
        sys.exit()

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    for i in range(len(clusters)): # pylint: disable=consider-using-enumerate
        if not os.path.exists(outdir + "/cluster_" + str(i)):
            os.mkdir(outdir + "/cluster_" + str(i))
        genes = list(clusters[i])
        with open(outdir + "/cluster_" + str(i) + "/genes_in_cluster",
                  "w", encoding = "utf8") as out:
            out.write("\n".join(genes))
        subset = edges[edges['Target'].isin(genes)]
        subset = subset[subset['Source'].isin(genes)]
        subset.to_csv(outdir + "/cluster_" + str(i) + "/cluster.csv",
                      index = False)

    #nx.draw(network, with_labels=True)
    #plt.show()




if __name__ == "__main__":
    main()
