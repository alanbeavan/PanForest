#!/usr/bin/env python3
"""Create and sql database for the network."""

import sys
import argparse
import sqlite3

def get_args():
    """Get user arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--edges", dest = "edges_file",
                        type = str, help = "Edges table")
    parser.add_argument("-n", "--nodes", dest = "nodes_file",
                        type = str, help = "Nodes table")
    parser.add_argument("-o", "--output", dest = "outfile",
                        type = str, help = "output file",
                        default = "network")
    args = parser.parse_args()
    if None in [args.edges_file, args.nodes_file]:
        parser.print_help(sys.stderr)
        sys.exit(0)
    return [args.edges_file, args.nodes_file, args.outfile]

def init_edges(infile, crsr, connection):
    """Initialise the Edges Table."""
    #statement for sql
    crsr.execute("DROP TABLE IF EXISTS edges")
    sql_command = """CREATE TABLE edges (
    EdgeID VARCHAR(255) PRIMARY KEY,
    Target_A VARCHAR(255) NOT NULL,
    Source_B VARCHAR(255) NOT NULL,
    InteractionType VARCHAR(2) NOT NULL,
    Weight DECIMAL(10,8) NOT NULL,
    P_B DECIMAL(10,8) NOT NULL,
    P_BGivenA DECIMAL(10,8) NOT NULL,
    P_BGivenNotA DECIMAL(10,8) NOT NULL,
    P_A DECIMAL(10,8) NOT NULL,
    P_AGivenB DECIMAL(10,8) NOT NULL,
    P_AGivenNotB DECIMAL(10,8) NOT NULL,
    P_AAndB DECIMAL(10,8) NOT NULL
    );"""
      
    # execute the statement
    crsr.execute(sql_command)

    #Read in table
    with open(infile) as in_file:
        lines = []
        for line in in_file:
            fields = line.rstrip("\n").split(",")
            values = ["--".join(fields[:2])]
            values.extend(fields)
            lines.append(values)

    #import statement
    connection. executemany("""
    INSERT INTO 
    edges 
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", lines[1:])

def init_nodes(infile, crsr, connection):
    """Initialise the Nodes Table."""
    #init
    crsr.execute("DROP TABLE IF EXISTS nodes")
    sql_command = """CREATE TABLE nodes (
    NodeID VARCHAR(255) PRIMARY KEY,
    Annotation VARCHAR(255),
    Description VARCHAR(255),
    Count INT,
    TruePositivesTest INT,
    FalsePositivesTest INT,
    FalseNegativesTest INT,
    TrueNegativesTest INT,
    ErrorRateTest DECIMAL(10,8),
    AccuracyTest DECIMAL(10,8),
    PrecisionClass1Test DECIMAL(10,8),
    PrecisionClass0Test DECIMAL(10,8),
    PrecisionAverageTest DECIMAL(10,8),
    RecallClass1Test DECIMAL(10,8),
    RecallClass0Test DECIMAL(10,8),
    RecallAverageTest DECIMAL(10,8),
    F1StatisticClass1Test DECIMAL(10,8),
    F1StatisticClass0Test DECIMAL(10,8),
    F1StatisticAverageTest DECIMAL(10,8),
    TruePositivesTrain INT,
    FalsePositivesTrain INT,
    FalseNegativesTrain INT,
    TrueNegativesTrain INT,
    ErrorRateTrain DECIMAL(10,8),
    AccuracyTrain DECIMAL(10,8),
    PrecisionClass1Train DECIMAL(10,8),
    PrecisionClass0Train DECIMAL(10,8),
    PrecisionAverageTrain DECIMAL(10,8),
    RecallClass1Train DECIMAL(10,8),
    RecallClass0Train DECIMAL(10,8),
    RecallAverageTrain DECIMAL(10,8),
    F1StatisticClass1Train DECIMAL(10,8),
    F1StatisticClass0Train DECIMAL(10,8),
    F1StatisticAverageTrain DECIMAL(10,8),
    DStatistic DECIMAL(10,8)
    );
    """
    #execute
    crsr.execute(sql_command)
    #Read in file
    with open(infile) as in_file:
        lines = []
        for line in in_file:
            fields = line.rstrip("\n").split(",")
            lines.append(fields)
    #import statement
    connection. executemany("""
    INSERT INTO 
    nodes
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", lines[1:])


def main():
    """Do the things."""
    edges_file, nodes_file, outfile = get_args()
  
    # connecting to the database
    connection = sqlite3.connect(outfile)
  
    # cursor
    crsr = connection.cursor()

    init_edges(edges_file, crsr, connection)
    init_nodes(nodes_file, crsr, connection)

    #crsr.execute("""select * from edges WHERE InteractionType = 'pp' LIMIT 10;""")
    #result = crsr.fetchall()
    #for x in result:
    #    print(x)
    


    # print statement will execute if there
    # are no errors
    print("Connected to the database")
  
    # close the connection
    connection.commit()
    connection.close()

if __name__ == "__main__":
    main()
