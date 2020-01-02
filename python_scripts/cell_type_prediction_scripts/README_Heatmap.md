#Cell Type Prediction 

This pipeline is used for generate heatmap to demonstrate peaks in multiple cell types. There are wig and bigwig approaches currently. 




##Required data
* peaks file: should contains peaks for segmenting significant/non-significant regions, it must represent each cell type peaks. Usually contains three columns: chromosome number (chr#), start position, end position. 
* SNPs file: a file with SNPs coordinates, should contains 6 columns: chromosome number (chr#), position of outside variant (start and end), outside variant, GWAS variant, significnat or not (1 represent for significant, 0 represent for non-significant). Example row: `chr2	230355955	230355956	rs10165012	rs10201872	1`
* cell type file: a file with one row only, which contains all cell types want to run in this pipeline. 
* wig/bigwig files: if use bigwig approach used, should give extra space, due to temporary bedgraph and wig files generated in the analysis. 


##Commands

###Wig approach
Example command for wig approach:

```
./PREDICT_celltype_peaksbypval_py_PV.sh \
-p /full/path/to/folder/ExampleData/Breast_and_Immune_collapsed_PVAL9_ext300_hg38.bed \
-s /full/path/to/folder/ExampleData/BINARY_Breast_Oct_hg38 \
-l test_cell_type_wig.txt \
-w /full/path/to/folder/ExampleData/WIGS_hg38/ \
-t wig \
-c 20 
-P /PubVer
```

###BigWig approach
The fllowing command used for bigwig file. 
For bigwig approach, [wigToBigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) should be pre-installed. 

```
./PREDICT_celltype_peaksbypval_py_PV.sh \
-p /full/path/to/folder/ExampleData/Breast_and_Immune_collapsed_PVAL9_ext300_hg38.bed \
-s /full/path/to/folder/ExampleData/BINARY_Breast_Oct_hg38 \
-l test_cell_type_bigwig.txt \
-w /full/path/to/folder/ExampleData/BW/ \
-t bigwig \
-c 20

```
 
Here is what each of the flag means:
####`-h/?`
These two flag use to get help information and exit process. Examle: `PREDICT_celltype_peaksbypval_py_PV.sh -h` or `PREDICT_celltype_peaksbypval_py_PV.sh -?`
####`-p`
This flag tells scripts which peaks file want to use and where can find it, must provide absolute path to peaks file. 
####`-s`
This flag gives which SNPs file will used in pipeline, absolute path to SNP file is required. 
####`-l`
This flag tels script where can find cell type list and the name of this file.
####`-d`
This flag is an optional argument, it tells script the distance from peak center want to used in this analysis, if it not provided 1500 will be used as default. 
####`-w`
This argument is required, in wig approach it tells script where the wig files is, or where can find bigwig file in bigwig approach. 
####`-t`
Requirement argument, tell script which approach will use in pipeline, wig approach or bigwig approach. 
####`-n`
An optional argument, tells script number of windows used in this analysis, if not given 40 will be used as default. 
####`-c`
An optional argument, if not given 12 will be used as default number. This flag tells script how many parallel jobs want to run same time, which is also same number of how many cores want to use. Maximum number of parallel jobs used in this pipeline is number of chromosomes, if the given number is greater than number of chromosomes want to run. Number of parallel jobs in this analysis is the cores number if cores number is smaller than number of chromosomes. 
####`-P`
This flag tell script the path to this pipeline (path to PREDICT_celltype_peaksbypval_py_PV.sh and Heatmap_Script_Files_PY folder), working directory will used as default if not provided. 

##Output files 

Directory:

####`HMW_peakFile_Distance`
This folder contains following subfolder:
#####`CellTypeName_OUTPUT`
The number of this format subfolder depends on how many cell types want to analysis. It contains subfolders `wig_Distance_NumWin` or  `bigwig_Distance_NumWin`, there are peaks files (separated based on chromosomes), searched wig records, z-score file. 
#####`Heatmap_Script_Files_PY`
A subdirectory contains all scripts required in pipeline, copied from the folder contains this pipeline scripts. 

`params_cellType_distance.txt` files generated in the analysis, each cell type has their own params file. 


##Other information
There are 22 chromosomes used in this pipeline, users can modify in the PREDICT_celltype_peaksbypval_py_PV.sh script to specify which chromosomes they want to use to test. 
