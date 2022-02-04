def preprocessDf(df, null_h, min_missing, min_present):
    """Modify the matrix so it's ready for random forest."""
    new_rows = []
    for row in df.index:
        new_row = []
        for field in row:
            new_row.append(str(field))
        new_rows.append(",".join(new_row))
    df.index = new_rows #convert row header to strings.
    df = df.fillna(0) #replace absent with 0
    df = df.replace({'.{2,}': '1'}, regex=True) #replace gene name with 1
    df = df.astype(int) #convert to numeric
    #remove genes with too much present
    df = df[df.sum(axis=1) < df.shape[1] - (min_missing -1)]
    #remove genes with too much absent
    df = df[df.sum(axis=1) > min_present - 1]
    if null_h:
        print("Randomising!!!!")
        for i in range(df.shape[0]):
            df.iloc[i] = random.sample(list(df.iloc[i]),
                                       len(list(df.iloc[i])))
    return df

def update_performance(table, i, y, y_train, y_test, y_trainP, y_testP):
    """Update the performance table."""
    #number of genomes present
    table['count'][i] = sum(list(y))
    cm_train = confusion_matrix(y_train, y_trainP)
    cm_test = confusion_matrix(y_test, y_testP)
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
    train_report = classification_report(y_train, y_trainP,
                                         output_dict = True,
                                         zero_division = 0)
    test_report = classification_report(y_test, y_testP,
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
    table['P1te'][i] = test_report['0']['precision']
    table['Pte'][i] = test_report['macro avg']['precision']
    #recalls - 1, 0 and average
    table['R1tr'][i] = train_report['1']['recall']
    table['R0tr'][i] = train_report['0']['recall']
    table['Rtr'][i] = train_report['macro avg']['recall']
    table['R1te'][i] = test_report['1']['recall']
    table['R1te'][i] = test_report['0']['recall']
    table['Rte'][i] = test_report['macro avg']['recall']
    #f1 scores - 1, 0 and average
    table['F1tr'][i] = train_report['1']['f1-score']
    table['F0tr'][i] = train_report['0']['f1-score']
    table['Ftr'][i] = train_report['macro avg']['f1-score']
    table['F1te'][i] = test_report['1']['f1-score']
    table['F0te'][i] = test_report['0']['f1-score']
    table['Fte'][i] = test_report['macro avg']['f1-score']


def fit_classifiers(df, imp, performance, ns, ng, ntrees, depth):
    """Fit a random forest classifier for all genes.

    Is it possible to vectorise this?!?!?!?!?!
    """
    df = df.transpose() #I think this is easier transposed
    for i in range(ng):
        print("gene number\t" + str(i+1) + "\tout of\t" + str(ng))
        y = df[df.columns[i]]
        x = df.drop([df.columns[i]], axis = 1)
        #now split the dataset into test and train - train with 75% in each
        #class
        X_train, X_test, y_train, y_test = train_test_split(x, y,
                                                            test_size=0.25,
                                                            stratify=y,
                                                            random_state=42)
        #random forest
        model = RandomForestClassifier(n_estimators = ntrees,
                                       max_depth = depth,
                                       max_features='sqrt',
                                       min_samples_split=2)
        train_fit = model.fit(X_train, y_train)
        y_pred_train = model.predict(X_train)
        y_pred_test = model.predict(X_test)
        #assess performance
        update_performance(performance, i, y, y_train, y_test,
                           y_pred_train, y_pred_test)
        #update importances (can I do tree and gene importances)
        imp[imp.columns[i]] = np.insert(model.feature_importances_, i, [0])
        if i % 1000 == 0 and i != 0:
            print("\nWriting importance and performance martices...\n")
            imp.to_csv("imp.csv")
            performance.to_csv("performance.csv")

