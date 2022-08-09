#!/usr/bin/env python3.6
"""Make an edge table with some extra stats."""

import argparse
import glob
import re
import sys
import pandas as pd
import rf_module as rf

def get_args():
    """Get user arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--performance", dest = "performance",
                        type = str, help = "Performance table")
    parser.add_argument("-d", "--dtable", dest = "d_table",
                        type = str, help = "D table")
    parser.add_argument("-o", "--output", dest = "output_file",
                        type = str, help = "Output table file")
    args = parser.parse_args()
    if None in [args.performance, args.d_table, args.output_file]:
        parser.print_help(sys.stderr)
        sys.exit(0)
    return [args.performance, args.d_table, args.output_file]


def main():
    """
    Write the nodes present in the graph with extra information.
    end result might be a table w/
    NodeID,Annotation,Wordy Annotation,[list of performance metrics].
    """
    #for now just do one cluster
    performance, d_table, output_file = get_args()
    d_stats = {}
    for line in rf.get_file_data(d_table)[1:]:
        fields = line.split("\t")
        d_stats[fields[0]] = fields[1]
    performance_lines = rf.get_file_data(performance)[1:]
    performance_lines = [line.replace("\"", "") for line in performance_lines]
    table_lines = ["NodeID,Annotation,Wordy Annotation,count,TPte,FPte,FNte,TNte,Ete,Ate,P1te,P0te,Pte,R1te,R0te,Rte,F1te,F0te,Fte,TPtr,FPtr,FNtr,TNtr,Etr,Atr,P1tr,P0tr,Ptr,R1tr,R0tr,Rtr,F1tr,F0tr,Ftr,D statistic"]
    for line in performance_lines:
        node_id = line.split(",")[0]
        try:
            table_lines.append(line + "," + d_stats[node_id])
            #print(line + "," + d_stats[node_id])
        except:
            #nodes not in any edges
            #print(node_id)
            continue

    with open(output_file, "w", encoding = "utf8") as out:
        out.write("\n".join(table_lines))


if __name__ == "__main__":
    main()
