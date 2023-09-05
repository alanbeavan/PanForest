# pangenome_rf
Pangenome analysis using random forests

**Author:** Alan Beavan

## Table of Contents
1. [Gather Materials](#1-gather-materials)
2. [Process Matrix](#2-process-matrix)
3. [Perform Random Forest](#3-perform-random-forest)
4. [Simplify the Importance Matrix (Optional)](#4-simplify-the-importance-matrix-optional)
5. [Convert to a Cytoscape/Gephi Format (Optional)](#5-convert-to-a-cytoscapegephi-format-optional)
6. [Direct the Edges of the Network (Optional)](#6-direct-the-edges-of-the-network-optional)
7. [Calculate the D Statistic for All Genes](#7-calculate-the-d-statistic-for-all-genes)
8. [Make an SQL Database (Optional)](#8-make-an-sql-database-optional)

---

## 1. Gather Materials
### a. Clone the GitHub Repository
Clone the GitHub repository using the following command:

git clone https://github.com/alanbeavan/pangenome_rf



### b. Infer the Pangenome
Infer the pangenome using panaroo. The file you need is `gene_presence_absence.csv`. You can find it in the GitHub repository. For detailed instructions, refer to their documentation.

### c. Install Python 3.9 or Above
Install Python 3.9 or a newer version. You can do this by installing from source locally or using Anaconda. To install using Anaconda, create a Python 3.10 environment:

conda create --name py10 python=3.10
conda activate py10


Then, install the required Python packages in your new Anaconda environment:

conda install pandas
conda install -c anaconda scikit-learn

If you installed a local version of Python, you can use `pip3` instead.

### d. Install Python Packages
Make sure you have installed the following Python packages in your environment:

- pandas
- scikit-learn

If you encounter any missing modules when running the program, install them using Anaconda.

## 2. Process Matrix
### a. Clean the Panaroo Presence-Absence Matrix
The Panaroo presence-absence matrix needs cleaning. Run `process_matrix.py` to clean it up:

python3 process_matrix.py <panaroo_presence_absence_matrix> <output_file>


My suggestion is to name the output file `collapsed_matrix.csv`. The program will generate six files as listed below:

- `collapsed_matrix.csv`: A cleaned version of the presence-absence matrix.
- `constant_genes.txt`: List of genes found in all genomes.
- `core_genes.txt`: List of genes found in >= 95% of genomes.
- `singletons.txt`: List of genes found in only one genome.
- `non_unique_genes.txt`: A table of gene family groupings.
- `non_unique_genomes.txt`: A table of genome groupings.

## 3. Perform Random Forest
### a. Familiarize with Options
First, familiarize yourself with the program options by running:

python3 pangenome_rf.py -h

This will display the program's usage.

### b. Run a Test
Run a short test to ensure everything is working:

python3 pangenome_rf.py -n 1 -d 1 -m collapsed_matrix.csv -pres 1 -abs 1 -o test


Check the output for progress information.

### c. Run the Random Forest Program
Run the random forest program with your chosen parameters. You need to decide on these main parameters:
- NTREES (number of trees in the forest)
- DEPTH (max depth of trees in the forest)
- MIN_PRESENT (minimum percentage of genomes featuring a gene)
- MIN_ABSENT (minimum percentage of genomes missing a gene)
- NTHREADS (number of parallel processes)

Example command:

python3 pangenome_rf.py -n <1000> -d <8> -m <collapsed_matrix.csv> -pres <1> -abs <1> -o <output>


Monitor the output directory for progress updates.

## 4. Simplify the Importance Matrix (Optional)
To simplify the importance matrix, set a threshold below which relationships should be ignored:

python3 simplify_imp.py <threshold> <input_imp.csv> <output_file>

The threshold should be chosen after reviewing the distribution of importances or network visualization.

## 5. Convert to a Cytoscape/Gephi Format (Optional)
To convert the matrix for use with Cytoscape and Gephi, run:

python3 convert_to_cytoscape.py <input imp.csv> <output_network.csv>

The output file will be a list of edges in the network.

## 6. Direct the Edges of the Network (Optional)
To indicate if genes are promoters or inhibitors, use the `direct_network.py` program:

python3 direct_network.py <input cytoscape csv> <collapsed_network.csv> <output>

## 7. Calculate the D Statistic for All Genes
To calculate the D statistic for genes, you'll need a rooted phylogeny and a list of genes. Run the following:

Rscript calculate_d.R -a <path> -t <treefile> -g <presence_absence_matrix> -c <ncores> -o <prefix>

The output is a file containing each gene and its D value.

## 8. Make an SQL Database (Optional)
To create an SQL database with statistics, follow these steps:

- Describe your nodes with performance metrics and D values:

python3 describe_nodes.py -p <PERFORMANCE> -d <D_TABLE> -o <OUTPUT_FILE>

- Describe your edges:

python3 describe_edges.py -m <MATRIX_FILE> -n <NETWORK_FILE> -o <OUTPUT_FILE>

- Create the SQL database:

python3 make_sql_database.py -e <EDGES_FILE> -n <NODES_FILE> -o <OUTFILE>

You now have an SQL database that you can query.

These instructions should help you with the Pangenome_rf.py workflow. Customize the commands and options based on your specific needs.
