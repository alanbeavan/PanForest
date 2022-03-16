#!/usr/bin/env python3.6
"""Convert importance matrices to Cyctoscape's preferred format."""

import sys
import pandas as pd

def get_args():
    """Get user arguments."""
    if len(sys.argv) != 3:
        print("USAGE: python3 convert_to_cytoscape.py infile outfile")
        sys.exit()
    return sys.argv[1:]

def main():
    """Columns: Source,Target,Interactiontype,Weight."""
    infile, outfile = get_args() # pylint: disable=unbalanced-tuple-unpacking
    lines = ["Source,Target,Interactiontype,Weight"]
    imp = pd.read_csv(infile, index_col = 0, header = 0)
    imp = imp.apply(pd.to_numeric)
    count = 0
    for index in imp.index:
        count += 1
        if count % 100 == 0:
            print(count)
        for col in imp[index].index:
            if imp[index][col] > 0:
                line = index.split(",")[0] + "," + col.split(",")[0] + ",pp," + str(imp[index][col])
                lines.append(line)
    with open(outfile, "w", encoding = "utf8") as out:
        out.write("\n".join(lines))


if __name__ == "__main__":
    main()
