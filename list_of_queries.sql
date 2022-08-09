/*
Here is a set of commands that can be run with
sqlite3 network ".read list_of_queries.sql"
Some will write to files - some are for analytical purposes only.

Alternatively, you can launch the sql interface using "sqlite3 network"
then running each query individually by copying and pasting.

Each query is prefaced with a comment describing it's purpose.
*/

--Get a list of genes that are predictable
.header ON
.mode csv
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/predictable_genes.txt
SELECT NodeID, Annotation, Description
FROM Nodes
WHERE ErrorRateTest <= 0.1
    AND F1StatisticClass1Test >= 0.9
    AND F1StatisticClass0Test >= 0.9
    AND DStatistic >= 0;

--lower threshold
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/lower_threshold_list.txt
SELECT NodeID, Annotation, Description
FROM Nodes
WHERE ErrorRateTest <= 0.2
    AND F1StatisticClass1Test >= 0.8
    AND F1StatisticClass0Test >= 0.8
    AND DStatistic >= 0;


--Get a list of genes that are predictable (not filtering for D)
.header ON
.mode csv
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/predictable_genes_ignore_D.txt
SELECT NodeID, Annotation, Description
FROM Nodes
WHERE ErrorRateTest <= 0.1
    AND F1StatisticClass1Test >= 0.9
    AND F1StatisticClass0Test >= 0.9;


--lower threshold with no D filtering
.header ON
.mode csv
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/lower_threshold_ignore_D.txt
SELECT NodeID, Annotation, Description
FROM Nodes
WHERE ErrorRateTest <= 0.2
    AND F1StatisticClass1Test >= 0.8
    AND F1StatisticClass0Test >= 0.8;
 

--Select all the edges of the network where the target is predictable and not phylogenetically correlated
.header ON
.mode csv
CREATE VIEW Temp_Edges 
    AS SELECT Target_A as Target,
        Source_B as Source,
        InteractionType,
        Weight,
        Description as Target_description
    FROM Edges E  INNER JOIN Nodes N ON E.Target_A = N.NodeID
    WHERE Target IN (
        SELECT NodeID
        FROM Nodes
        WHERE ErrorRateTest <= 0.1
            AND F1StatisticClass1Test >= 0.9
            AND F1StatisticClass0Test >= 0.9
            AND DStatistic > 0);
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/conservative_network.csv
SELECT Target,
    Source,
    InteractionType,
    Weight,
    Target_description,
    Description as Source_description
FROM Temp_Edges E INNER JOIN Nodes N ON E.Source = N.NodeID;




/*SELECT Target_A as Target,
    Nodes.Description where Nodes.NodeID = Target_A as Target_description
    Source_B as Source,
    Nodes.Description where Nodes.NodeID = Source_B as Source_description
    InteractionType,
    Weight
FROM Edges
WHERE Target IN (
    SELECT NodeID
    FROM Nodes
    WHERE ErrorRateTest <= 0.1
        AND F1StatisticClass1Test >= 0.9
        AND F1StatisticClass0Test >= 0.9
        AND DStatistic > 0);
*/
--obtain the commensal interactions
.header ON
.mode csv
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/commensal_edges.csv
SELECT * FROM edges WHERE Target_A IN (
    SELECT NodeID FROM Nodes WHERE ErrorRateTest <= 0.1 AND F1StatisticClass1Test >= 0.9 AND F1StatisticClass0Test >= 0.9 AND DStatistic > 0)
AND ((P_AGivenB >= 0.95 AND P_BGivenA <= 0.95 AND P_AGivenNotB >= 0.05 AND P_AGivenNotB >= P_A/2)
OR (P_BGivenA >= 0.95 AND P_AGivenB <= 0.95 AND P_BGivenNotA >= 0.05 AND P_BGivenNotA >= P_B/2))
AND InteractionType = "pp";

--More liberal commensals
.header ON
.mode csv
CREATE VIEW Commensals
    AS SELECT * FROM edges WHERE Target_A IN (
        SELECT NodeID FROM Nodes WHERE ErrorRateTest <= 0.1 AND F1StatisticClass1Test >= 0.9 AND F1StatisticClass0Test >= 0.9 AND DStatistic > 0)
    AND ((P_AGivenB >= 0.99 AND P_BGivenA <= 0.99 AND P_AGivenNotB >= 0.01 AND P_AGivenNotB >= P_A/5)
    OR (P_BGivenA >= 0.99 AND P_AGivenB <= 0.99 AND P_BGivenNotA >= 0.01 AND P_BGivenNotA >= P_B/5))
    AND InteractionType = "pp";


.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/liberal_commensal_edges.csv
SELECT * from Commensals;


--obtain the mutualistic interactions
.header ON
.mode csv
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/mutual_edges.csv
SELECT * FROM edges
WHERE Target_A IN (
    SELECT NodeID
    FROM Nodes
    WHERE ErrorRateTest <= 0.1
    AND F1StatisticClass1Test >= 0.9
    AND F1StatisticClass0Test >= 0.9
    AND DStatistic > 0)
AND edges.EdgeID in (
    SELECT Source_B || '--' || Target_A
    FROM edges
    WHERE Target_A IN (
        SELECT NodeID
        FROM Nodes
        WHERE ErrorRateTest <= 0.1
        AND F1StatisticClass1Test >= 0.9
        AND F1StatisticClass0Test >= 0.9
        AND DStatistic > 0))
AND EdgeID not in (SELECT EdgeID FROM Commensals)
AND Target_A || '--' || Source_B NOT IN (SELECT EdgeID FROM Commensals);

--Obtain the competition reactions
.header ON
.mode csv
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/competition_edges.csv
SELECT * FROM edges
WHERE Target_A IN (
    SELECT NodeID
    FROM Nodes
    WHERE ErrorRateTest <= 0.1
    AND F1StatisticClass1Test >= 0.9
    AND F1StatisticClass0Test >= 0.9
    AND DStatistic > 0)
AND InteractionType = "nn";



--Obtain the most connected nodes (Source then Target, +ve then -ve)
.header ON
CREATE VIEW MostTargettedPos
AS SELECT Source_B AS Source,
        P_A AS Freq,
        COUNT(*) AS Edges_emerging
    FROM Edges
    WHERE InteractionType = 'pp' 
    AND Target_A IN (
        SELECT NodeID
        FROM Nodes
        WHERE ErrorRateTest <= 0.1
        AND F1StatisticClass1Test >= 0.9
        AND F1StatisticClass0Test >= 0.9
        AND DStatistic > 0)
    GROUP BY Source_B
    ORDER BY Edges_emerging DESC
limit 300;

.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/most_predicting_positive.csv
SELECT NodeID,
    Annotation,
    Description,
    Freq,
    Edges_emerging
FROM MostTargettedPos M INNER JOIN Nodes N ON M.Source = N.NodeID;


.header ON
CREATE VIEW MostTargetingPos AS
SELECT Target_A AS Target,
        P_A AS Freq,
        COUNT(*) AS Edges_targeting
    FROM Edges
    WHERE InteractionType = 'pp' 
    AND Target_A IN (
        SELECT NodeID
        FROM Nodes
        WHERE ErrorRateTest <= 0.1
        AND F1StatisticClass1Test >= 0.9
        AND F1StatisticClass0Test >= 0.9
        AND DStatistic > 0)
    GROUP BY Target_A
    ORDER BY Edges_targeting DESC
    limit 300;
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/most_predicted_positive.csv
SELECT NodeID,
    Annotation,
    Description,
    Freq,
    Edges_targeting
FROM MostTargetingPos M INNER JOIN Nodes N ON M.Target = N.NodeID;


.header ON
CREATE VIEW MostTargettedNeg AS
    SELECT Source_B AS Source,
        P_A AS Freq,
        COUNT(*) AS Edges_emerging
    FROM Edges
    WHERE InteractionType = 'nn' 
    AND Target_A IN (
        SELECT NodeID
        FROM Nodes
        WHERE ErrorRateTest <= 0.1
        AND F1StatisticClass1Test >= 0.9
        AND F1StatisticClass0Test >= 0.9
        AND DStatistic > 0)
    GROUP BY Source_B
    ORDER BY Edges_emerging DESC
    limit 300;
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/most_predicting_negative.csv

SELECT NodeID,
    Annotation,
    Description,
    Freq,
    Edges_emerging
FROM MostTargettedNeg M INNER JOIN Nodes N ON M.Source = N.NodeID;


.header ON
CREATE VIEW MostTargetingNeg AS
    SELECT Target_A AS Target,
        P_A AS Freq,
        COUNT(*) AS Edges_targeting
    FROM Edges
    WHERE InteractionType = 'nn' 
    AND Target_A IN (
        SELECT NodeID
        FROM Nodes
        WHERE ErrorRateTest <= 0.1
        AND F1StatisticClass1Test >= 0.9
        AND F1StatisticClass0Test >= 0.9
        AND DStatistic > 0)
    GROUP BY Target_A
    ORDER BY Edges_targeting DESC
    limit 300;
.once /Users/ab17362/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Post_doc_notts/sqlite_practice/most_predicted_negative.csv

SELECT NodeID,
    Annotation,
    Description,
    Freq,
    Edges_targeting
FROM MostTargetingNeg M INNER JOIN Nodes N ON M.Target = N.NodeID;


