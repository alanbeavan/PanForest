"""
rf_module.py

Module used for the analysis of pangenomes using random forests.
"""
import random
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report

def get_file_data(filename):
    """Stores the lines of the program with name filename as a list."""
    with open(filename, encoding = "utf8") as in_file:
        lines = []
        for line in in_file:
            lines.append(line.rstrip("\n"))
    return lines

def init_tables(table):
    """Initialise the performance and importance tables."""
    n_g = table.shape[0] #number of genes (gene families)
    imp = pd.DataFrame(0.0, index = np.arange(n_g),
                       columns = table.index.values)
    imp.index = table.index.values
    metrics = ['count', 'TPte', 'FPte', 'FNte', 'TNte', 'Ete', 'Ate', 'P1te',
               'P0te', 'Pte', 'R1te', 'R0te', 'Rte', 'F1te', 'F0te', 'Fte',
               'TPtr','FPtr', 'FNtr', 'TNtr', 'Etr', 'Atr', 'P1tr', 'P0tr',
               'Ptr', 'R1tr', 'R0tr', 'Rtr', 'F1tr', 'F0tr', 'Ftr']
    performance = pd.DataFrame(0.0, index = np.arange(n_g), columns = metrics)
    performance.index = table.index.values
    return imp, performance

def preprocess_df(table, null_h, min_missing, min_present):
    """Modify the matrix so it's ready for random forest."""
    new_rows = []
    for row in table.index:
        new_row = []
        for field in row:
            new_row.append(str(field))
        new_rows.append(",".join(new_row))
    table.index = new_rows #convert row header to strings.
    table = table.fillna(0) #replace absent with 0
    table = table.replace({'.{2,}': '1'}, regex=True) #convert to 1
    table = table.astype(int) #convert to numeric
    #remove genes with too much present
    table = table[table.sum(axis=1) < table.shape[1] - (min_missing -1)]
    #remove genes with too much absent
    table = table[table.sum(axis=1) > min_present - 1]
    if null_h:
        print("Randomising!!!!")
        for i in range(table.shape[0]):
            table.iloc[i] = random.sample(list(table.iloc[i]),
                                       len(list(table.iloc[i])))
    return table

def update_performance(table, i, y_sets):
    """Update the performance table."""
    #number of genomes present
    table['count'][i] = sum(list(y_sets[0]))
    cm_train = confusion_matrix(y_sets[1], y_sets[3])
    cm_test = confusion_matrix(y_sets[2], y_sets[4])
    #stats for training set
    table['TPtr'][i] = cm_train[1,1]
    table['FPtr'][i] = cm_train[0,1]
    table['TNtr'][i] = cm_train[0,0]
    table['FNtr'][i] = cm_train[1,0]
    #stats for test set
    table['TPte'][i] = cm_test[1,1]
    table['FPte'][i] = cm_test[0,1]
    table['TNte'][i] = cm_test[0,0]
    table['FNte'][i] = cm_test[1,0]
    #reports for recall, precision, f1, accuracy
    train_report = classification_report(y_sets[1], y_sets[3],
                                         output_dict = True,
                                         zero_division = 0)
    test_report = classification_report(y_sets[2], y_sets[4],
                                        output_dict = True,
                                        zero_division = 0)
    table['Etr'][i] = 1 - train_report['accuracy']
    table['Ete'][i] = 1 - test_report['accuracy']
    table['Atr'][i] = train_report['accuracy']
    table['Ate'][i] = test_report['accuracy']
    #precisions - 1, 0 and average
    table['P1tr'][i] = train_report['1']['precision']
    table['P0tr'][i] = train_report['0']['precision']
    table['Ptr'][i] = train_report['macro avg']['precision']
    table['P1te'][i] = test_report['1']['precision']
    table['P0te'][i] = test_report['0']['precision']
    table['Pte'][i] = test_report['macro avg']['precision']
    #recalls - 1, 0 and average
    table['R1tr'][i] = train_report['1']['recall']
    table['R0tr'][i] = train_report['0']['recall']
    table['Rtr'][i] = train_report['macro avg']['recall']
    table['R1te'][i] = test_report['1']['recall']
    table['R0te'][i] = test_report['0']['recall']
    table['Rte'][i] = test_report['macro avg']['recall']
    #f1 scores - 1, 0 and average
    table['F1tr'][i] = train_report['1']['f1-score']
    table['F0tr'][i] = train_report['0']['f1-score']
    table['Ftr'][i] = train_report['macro avg']['f1-score']
    table['F1te'][i] = test_report['1']['f1-score']
    table['F0te'][i] = test_report['0']['f1-score']
    table['Fte'][i] = test_report['macro avg']['f1-score']


def fit_classifiers(table, results, params, output, checkpoint):
    """
    Fit a random forest classifier for all genes.
    results = [imp, performance]
    params = [ntrees, depth, nthreads]
    """
    n_g = table.shape[0]
    table = table.transpose() #I think this is easier transposed
    start = 0
    if checkpoint != 0:
        results[0] = pd.read_csv(output +"/imp.csv", header = 0, index_col = 0)
        results[1] = pd.read_csv(output + "/performance.csv", header = 0, index_col = 0)
        start = checkpoint

    for i in range(start, n_g):
        print("gene number\t" + str(i+1) + "\tout of\t" + str(n_g))
        y_all = table[table.columns[i]]
        x_all = table.drop([table.columns[i]], axis = 1)
        #now split the dataset into test and train - train with 75% in each
        #class
        # datasets = x_train, x_test, y_train, y_test
        datasets = train_test_split(x_all, y_all,
                                    test_size=0.25, stratify=y_all)
        #random forest
        model = RandomForestClassifier(n_estimators = params[0],
                                       max_depth = params[1],
                                       max_features='sqrt',
                                       min_samples_split=2,
                                       n_jobs=params[2])
        model.fit(datasets[0], datasets[2])
        y_pred_train = model.predict(datasets[0])
        y_pred_test = model.predict(datasets[1])
        #assess performance
        update_performance(results[1], i,
                           [y_all, datasets[2], datasets[3],
                            y_pred_train, y_pred_test])
        #update importances
        results[0][results[0].columns[i]] = \
                np.insert(model.feature_importances_, i, [0])

        if i % 1000 == 0 and i != 0:
            print("\nWriting importance and performance martices...\n")
            results[0].to_csv(output + "/imp.csv")
            results[1].to_csv(output + "/performance.csv")
    return results
