# tSpace algorithm for unsupervised deteremination of multiple developmental branches in single cell data

By Denis Dermadi

Description

tSpace is an algorithm for trajectory inference described in the publication (https://doi.org/10.1101/336313).

Originally, it was developed for single cell analysis, however it can be applied on any type of large data.

Here we focus on FACS, CyTOF or single cell (sc) RNAseq data. MATLAB version of tSpace requires previously cleaned data of all the noise/artifacts. Data does not need to be previously transformed. tSpace script provides asinh and logicle transformations for FACS and CyTOF, and log2+1 transformation and scaling for scRNAseq data. Expression matrix contains cells in the rows, and measured proteins/genes in the columns.

For the first analysis we recommend to use default values, with the exception of waypoints (wp) which for smaller complexity data can be reduced to e.g. here 15. Usually 20 is a good number. MATALB version of tSpace supports for visualization PCA. However, users can visualize trajectory matrix using UMAP in e.g. Python or R.

## Commentary on tSpace parameters

### Waypoints and number of trajectories

In the accompanying publication, we examine (i) the effects of varying T (trajectory number) and effect of waypoints (ii) and the effects of graph number (G), K (L is maintained at 0.75K as default) and metric parameters. We demonstrate the power of waypoints (WP) and the importance of the number of trajectories on tSpace performance. WP are crucial element for tSpace to perform well. Based on our experience, the complexity of the dataset dictates number of WP; so far, we have used anything from 10 to 30 waypoints. 100-200 trajectories are needed to reveal the details of branching, however, for final analysis we suggest to calculate preferably larger number of trajectories (~1000) or ground truth if working with fewer than 10000 cells.

### Comments on K and L

Furthermore, we compare the effect of the number of sub-graphs, and varying K, L and different metrics (Euclidean and Pearson) and demonstrate the robustness of the output to K, the number of neighbors in the KNN graph, as long as K is kept reasonably low (but high enough for connectivity of the KNN graph). Intuitively large K, by increasing the number of ‘paths’ between cells, can lead to unwanted connections (short circuits) to developmentally more distant or altogether unrelated cells. We suggest to keep K small, and set the default to K=20. The parameter L defines the subset of K connections around each cell to be preserved in each of the sub-graphs, thus L must be < K. Keeping K and L small increases the sparseness of the graphs (reducing aberrant paths); but the ratio of L to K determines the independence of sub-graphs (which is also important to reduce the contribution of short-circuits). If L is kept as NULL tSpace will determine it as 0.75`*`K, a value that works well in all datasets we have analyzed. User can override that default by simply specifying L value.

### Comments on number of graphs

Increased number of sub-graphs (graph parameter) can potentially contribute to noisy connections between cell types that are not developmentally related, however from developmental point of view, irrespective of choice of graph numbers cells are ordered correctly.

### Comments on metrics

In our analysis we used Euclidean and Pearson correlation (correlation) metric, all commonly used in single cell analysis. Selection of metrics is dependent on the data type. For example, Euclidean distances may over-emphasize the contribution of phenotypic markers that are very highly expressed, unless markers are scaled to a similar range prior to tSpace. Pearson metric intrinsically compensates for this, focusing on the profile shapes rather than magnitudes.


# Installation

Copy src directory and tspace.m file into a tspace directory within your MATALB working directory and add to path directory and subdirectories. 

**Note**: 
tSpace depends on 
Statistics and Machine Learning Toolbox, 
Parallel Computing Toolbox, 
Bioinformatics Toolbox


# How to use the algorithm

Type in command window
`run tspace.m`

and follow the instructions. 
Function will export csv files in the working directory, which can be visualized and processed further, see more [detailed description](http://denisdermadi.com/tspace-trajectory-inference-algorithm).

# tcell_demo

A csv file containing logicle transformed FACS data of 10000 T cells. This file is mean't for testing MATLAB function.
