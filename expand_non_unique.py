#!/usr/bin/env python3.6
"""
Expand the non-unique gene presesence absence patterns so that each family
group is replaced with the genes that belong to it. The resulting network
will no longer feature edge weights because new connections between genes
of identical PA pattern will have weights that cannot be defined in the
same way.

Args:
    network:        The network to be expanded
    groups:         The set of non-unique gene PA patterns outputted by 
                    process_matrix
    include_new:    (optional) An indication that the new network should 
                    include family group expansions that do not currently
                    feature if the supplied network. i.e. if a family group
                    featuring 2 genes had no edges in the suplied network, a
                    new connected component would appear in the output
                    network featuring these 2 genes as co-occuring but with
                    no other connections.
    node_table:     (optional) If include_new = True, the full node table
                    outputted by describe_nodes.py must be passed to this
                    program so that family groups not currently in the network
                    can be assessed for whether they pass the D score
                    threshold. Only genes present in the node_table will be
                    added even with -i so make sure it has everything you need
    D_cutoff:       (optional) If include_new = True, the D score minimum for
                    new family groups to be included in the output network.

"""

import argparse
import sys
import rf_module as rf

def get_args():
    """Get settings from user arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--network-in", type = str,
                        help = "File name of network to modify",
                        dest = "net_file")
    parser.add_argument("-g", "--groups", type = str,
                        help = "Non-unique gene families file as output by \
                                process_matrix.py",
                        dest = "groups_file")
    parser.add_argument("-i", "--include-new",
                        help = "include links between family groups not curren\
                                tly in the network",
                        dest = "include", action = "store_true")
    parser.add_argument("-t", "--node-table", type = str,
                        help = "The full node table output by describe_nodes.p\
                                y if -i",
                        dest = "node_table")
    parser.add_argument("-d", "--d-cutoff", type = float,
                        help = "minimum value of D for which a gene can be\
                                 included in the new network if -i",
                        default = 0.0, dest = "d_min")
    parser.add_argument("-o", "--out", type = str,
                        help = "output network file name",
                        default = "expanded_network.csv", dest = "outfile")
    args = parser.parse_args()
    if None in [args.net_file, args.groups_file]:
        parser.print_help(sys.stderr)
        sys.exit(0)
    return [args.net_file, args.groups_file, args.include, args.node_table,
            args.d_min, args.outfile]

def self_match(g_list):
    """return the combinations of genes pairs genes in the same group."""
    matches = []
    for gene in g_list:
        for gene1 in g_list:
            if gene != gene1:
                matches.append(gene + "," + gene1 + ",pp")
    return matches

def main():
    """Expand family groups and write new output."""
    network, groups, include, node_table, d_min, outfile = get_args()
    new_net = []
    groups_dict = {}
    for line in rf.get_file_data(groups):
        fields = line.split("\t")
        genes = []
        gene_fields = fields[1].split(",")
        i = 0
        while i < len(gene_fields):
            genes.append(gene_fields[i])
            i += 3
        groups_dict[fields[0]] = genes
    
    done_fams = []
    for line in rf.get_file_data(network)[1:]:
        if "family" in line:
            fields = line.split(",")
            int_type = fields[2]
            if "family" in fields[0] and "family" in fields[1]:
                for tgene in groups_dict[fields[0]]:
                    for sgene in groups_dict[fields[1]]:
                        new_net.append(tgene + "," + sgene + "," + int_type)
                new_net.extend(self_match(groups_dict[fields[0]]))
                new_net.extend(self_match(groups_dict[fields[1]]))
                if include:
                    if fields[0] not in done_fams:
                        done_fams.append(fields[0])
                    if fields[1] not in done_fams:
                        done_fams.append(fields[1])
            elif "family" in fields[0]:
                for tgene in groups_dict[fields[0]]:
                     new_net.append(tgene + "," + fields[1] + "," + int_type)
                new_net.extend(self_match(groups_dict[fields[0]]))
                if include and fields[0] not in done_fams:
                    done_fams.append(fields[0])

            elif "family" in fields[1]:
                for sgene in groups_dict[fields[1]]:
                    new_net.append(fields[0] + "," + sgene + "," + int_type)
                new_net.extend(self_match(groups_dict[fields[1]]))
                if include and fields[1] not in done_fams:
                    done_fams.append(fields[1])
        else:
            new_net.append(",".join(line.split(",")[:3]))
    
    if include:
        for line in rf.get_file_data(node_table):
            if "family" in line:
                fields = line.split(",")
                if "family" in fields[0]  and fields[0] not in done_fams:
                    if float(fields[-1]) > d_min:
                        new_net.extend(self_match(groups_dict[fields[0]]))
    with open(outfile, "w", encoding = "utf8") as out:
        out.write("Target,Source,InteractionType\n")
        out.write("\n".join(new_net) + "\n")
                    

if __name__ == "__main__":
    main()
