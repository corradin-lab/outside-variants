#!/bin/bash
#FileName: PREDICT_celltype_peaksbypval_py_PV.sh
#Modified by: Yanwei Song 
#Purpose: 

#Input:

### IMPORTANT NOTES:
#### NOTE 1: if HMW_$peaks is not unique files will be overwritten, run from unique folder for each analysis is safer 
#### NOTE 2: uses 4 processors for HMW so use with bsub -n 4
#### NOTE 3: runs heatmap linearly
#### NOTE 4: peaks and snps need MATCHING coordinate systems hg38 for both
###### NOTE 5: DON't run wilcox and ks on same list in parallel from same directory files will get overwritten, run KS than run wilcox with NO make heatmap

######	Run heatmapwindow yay or nay

####### Peaks list
####### binary snp file coordinates
####### path to directory with cell type wigs
####### file list of cell type names while read a  done < cell type list

## peaks is peaks for segmenting sig/nsig regions: IMPORTANT MUST REPRESENT EACH CELL TYPE PEAKS, DUE TO ZSCORE IE IF NEURONAL SAMPLE RUN WITH ONLY BLOOD PEAKS INFLATE NEURONAL SIGNAL
## $snps is binary file with snps coordinates example row: chr2	230355955	230355956	rs10165012	rs10201872	1
## $cell_types is file with one row only, input for for loop, listing sample names exactly
 
## $run_heatmap --> binary yes if calculate HMW script 
## $dis is distance from peak center, ie 2500 -> +/-2500 5kb window for hmw script
## $path is path to WIG files for HMW script : example : /lab/corradin_data/OutVar_by_celltype/IMPUTE_WTCCC/2018_outvars/Nov18_Beyond/HG38_nov18/wigs/
## $r is ksonesided OR wilcox to call either of these scripts in /lab/corradin_data Aggregate_ksonesided_paired_greater_tests.r or Aggregate_wilcox_paired_greater_tests.r
## $Ronly if last variable is Ronly (exactly) it will skip all things EXCEPT running the Rtests, good for when already ran ks and now want wilcox
#example:  /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/PubVer/PREDICT_celltype_peaksbypval_py_PV.sh -p /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/Breast_and_Immune_collapsed_PVAL9_ext300_hg38.bed -s /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/BINARY_Breast_Oct_hg38 -l /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/test_cell_type.txt -w /lab/corradin_data/permanent/WIGS_COMBINED/hg38/ -t wig -c 20 -P /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/PubVer

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-hpsldwtncP] [-p PeaksFile] [FILE]...
Hello Users! 
This script is used to predict peaks of multiple cell types based on p-value pipeline. 

### IMPORTANT NOTES:
#### NOTE 1: if HMW_$peaks is not unique files will be overwritten, run from unique folder for each analysis is safer 
#### NOTE 2: Max processors for HMW can use is max number of chromosomes number (e.g. 22 for humans)
#### NOTE 3: runs heatmap linearly
#### NOTE 4: peaks and snps need MATCHING coordinate systems hg38 for both
#### NOTE 5: DON't run wilcox and ks on same list in parallel from same directory files will get overwritten, run KS than run wilcox with NO make heatmap
#### NOTE 6: Must give 5 arguments: peaks file name, SNPs file name, cell type list, path to wig/bigwig file(s), file type
#### NOTE 7: For bigwig approach, user should to be noticed this pipeline will generate some extra files (like .bedgraph and .wig files), so it requires extra spaces when scripts running. 

###### Run heatmapwindow yay or nay
###### Peaks list
###### binary snp file coordinates
###### path to directory with cell type wigs
 
###### Example(submited on server): bsub "/lab/corradin_data/permanent/CELL_TYPE_PREDICTION/PubVer/PREDICT_celltype_peaksbypval_py_PV.sh -p /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/Breast_and_Immune_collapsed_PVAL9_ext300_hg38.bed -s /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/BINARY_Breast_Oct_hg38 -l /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/test_cell_type.txt -w /lab/corradin_data/permanent/WIGS_COMBINED/hg38/ -t wig -c 20 -P /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/PubVer"


    -h/? 			Display this help and exit.
    -p PeaksFile	File name. Peaks for segmenting sig/nsig regions: IMPORTANT MUST REPRESENT EACH CELL TYPE PEAKS, DUE TO ZSCORE IE IF NEURONAL SAMPLE RUN WITH ONLY BLOOD PEAKS INFLATE NEURONAL SIGNAL. (whole path to peaks file)
    -s SNPs		File name. Is binary file with snps coordinates example row: chr2	230355955	230355956	rs10165012	rs10201872	1
    -l CellTypeList	File name. Is a file with one row only, contains all cell types want to use in this pipeline.
    -d distance		Distance from peak center, ie 2500 -> +/-2500 5kb window for hmw script. Optional argument, if no value given, 1500 will be used as default.    		
    -w /../../		Whole path to wig (or bigwig) files. 
    -t FileType		File type user want to use, wig or bigwig. 
    -n NumberWindow	Number of windows, optional argumet, if not given 40 will be used as default. 
    -c Cores        The number of subprocessors user want to use, related to the number of cores are available to use, default is 12 
    -P PATH         Path to where user put this script and relative folder(Heatmap_Script_Files_PY), if not provided will use current working directory. 

EOF
1>&2
exit 1
}

while getopts ":hp:s:l:d:w:t:n:c:P:" opt
do
    case "${opt}" in
        h|\?)
            show_help
            exit 0
            ;;
        p)  peaksP=${OPTARG}
            ;;
        s)  snpsP=${OPTARG}
            ;;
        l)  cell_list=${OPTARG}
            ;;
        d)  dis=${OPTARG}
            ;;
        w)  wig_path=${OPTARG}
            ;;
        t)  type=${OPTARG}
            ;;
        n)  num_windows=${OPTARG}
            ;;
        c)  num_cores=${OPTARG}
            ;;
        P)  path=${OPTARG}
            ;;
        
    esac
done
shift $((OPTIND-1))   # Discard the options and sentinel 


# if there no argument give, provide usage information
if [ $OPTIND -eq 1 ]
# if [ $OPTIND -lt 19 ]
then
    show_help
    exit 0
fi
    
# check if required argument provided:
if [[ -z $peaksP || -z $snpsP || -z $cell_list || -z $wig_path ||  -z $type ]]
then
    echo "
    ERROR:
    Please check if peaks, snps, cell list, path to wig(or bigwig) files or file type (bigwig/wig) is provided!!
    
    "
    show_help
    exit 0
fi

# if number of windows not given, use 40 as default
if [ -z "$num_windows" ]
then
	num_windows=40
fi
# if distance from peak center not given, use 1500 as default
if [ -z "$dis" ]
then
	dis=1500
fi
# if number of cores want to use not given, use 12 as default
if [ -z "$num_cores" ]
then
	num_cores=12
fi
# if path to Heatmap_Script_Files_PY not given, use current directory as default
if [ -z "$path" ]
then
	path=$(pwd)
fi

#a list extract from cell list, input for for loop, listing sample names exactly
cell_types=$(<$cell_list)

#convert all letters to lowercase 
# r=`echo "$r" | tr '[:upper:]' '[:lower:]'`


#--------------Main pipeline start from here--------------
# if [ $Ronly != "Ronly" ]
# then

echo "Run pipeline"

# if [ $run_heatmap = "yes" ]
# then

#split peaks and snps path just keep file name, in case peaks and snps files are not in the working directory
peaks=`echo ${peaksP##*/}`
snps=`echo ${snpsP##*/}`

## Prep for heatmap windows
mkdir HMW_"$peaks"_"$dis"
# cp /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/params_test.txt HMW_"$peaks"_"$dis"/   #########why? really necessary?
# cp -r /lab/corradin_data/permanent/CELL_TYPE_PREDICTION/PubVer/Heatmap_Script_Files_PY/ HMW_"$peaks"_"$dis"/
#get path to shell scrip
if [ ${path: -1} == "/" ]
then
cp -r "$path"Heatmap_Script_Files_PY/ HMW_"$peaks"_"$dis"/
else
cp -r "$path"/Heatmap_Script_Files_PY/ HMW_"$peaks"_"$dis"/
fi
# cp $peaks HMW_"$peaks"_"$dis"/    ############may can remove after analysis
# cp $peaksP HMW_"$peaks"_"$dis"/    ############may can remove after analysis

## prep params file for hmw
for cells in $cell_types
do
	echo "chromosomes:1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22" > HMW_"$peaks"_"$dis"/params_"$cells"_"$dis".txt
	# echo "absolute_path_to_peak_file:$peaks" >> HMW_"$peaks"_"$dis"/params_"$cells"_"$dis".txt
	echo "absolute_path_to_peak_file:$peaksP" >> HMW_"$peaks"_"$dis"/params_"$cells"_"$dis".txt
	echo "absolute_path_to_sample_dirs:$wig_path" >> HMW_"$peaks"_"$dis"/params_"$cells"_"$dis".txt
	echo "sample_dirs:$cells" >> HMW_"$peaks"_"$dis"/params_"$cells"_"$dis".txt
	echo "calculate_peak_center?:yes" >> HMW_"$peaks"_"$dis"/params_"$cells"_"$dis".txt
	echo "num_windows:$num_windows" >> HMW_"$peaks"_"$dis"/params_"$cells"_"$dis".txt
	echo "dist_from_peak_in_each_direction:$dis" >> HMW_"$peaks"_"$dis"/params_"$cells"_"$dis".txt
	echo "file_type:$type" >> HMW_"$peaks"_"$dis"/params_"$cells"_"$dis".txt #change here if you want to use bigwig files instead of wig files 
done 
cur_dir=$(pwd)    
cd HMW_"$peaks"_"$dis"

#get working directory
dir=$(pwd)
# echo "$dir"/Heatmap_Script_Files_PY/runHeatmapWindowMedians_YS.py
for cells in $cell_types
do	    
# perl runHeatmapWindowMedians.pl params_"$cells"_"$dis".txt 4
python "$dir"/Heatmap_Script_Files_PY/runHeatmapWindowMedians_YS.py "$dir"/params_"$cells"_"$dis".txt "$num_cores"
done


echo "Heatmap windows complete"

cd ../

# else
# 	echo "heatmap script not run"
	
# fi


####### make peaks binary based on SNPs

### peaks to midpoint command: 

cat "$peaksP" | awk '{printf "%s\t%d\t%d\t%s\t%0.1f\n", $1,$2,$3,$1,$2+($3-$2)/2}' | sed 's/\.0//' | sed 's/\t/_/4' > Mdpts_"$peaks"

### BINARY PEAKS T2 based on matched pvalue thresholds

intersectBed -a "$snpsP" -b Mdpts_"$peaks" -wa -wb | awk '($6==1)' | awk '{OFS="\t"}{print $7,$8,$9,$10,$5}' | sort -u > SIG_PEAKS_"$peaks"_"$snps"
intersectBed -a "$snpsP" -b Mdpts_"$peaks" -wa -wb | awk '{OFS="\t"}{print $7,$8,$9,$10,$5}' | sort -u  > TESTED_PEAKS_"$peaks"_"$snps"
cat SIG_PEAKS_"$peaks"_"$snps" TESTED_PEAKS_"$peaks"_"$snps"| sort -d | uniq -c | awk '($1==2)' | awk '{OFS="\t"}{print $2,$3,$4,$5,$6,$7,1}' > BINARY_PEAKS_"$peaks"_"$snps"
cat SIG_PEAKS_"$peaks"_"$snps" TESTED_PEAKS_"$peaks"_"$snps" | sort -d | uniq -c | awk '($1==1)' | awk '{OFS="\t"}{print $2,$3,$4,$5,$6,$7,0}' >> BINARY_PEAKS_"$peaks"_"$snps"


####### GET AVERAGES for sig and nsig peaks 

## CONNECT BINARY FILES TO zscores

wc -l BINARY_PEAKS_"$peaks"_"$snps"

cat BINARY_PEAKS_"$peaks"_"$snps" | awk '($6==1)' | cut -f 5 | sort -u > BINsnps_"$peaks"_"$snps"
cat BINARY_PEAKS_"$peaks"_"$snps" | awk '($6==0)' | cut -f 5 | sort -u >> BINsnps_"$peaks"_"$snps"
cat BINsnps_"$peaks"_"$snps" | sort -d | uniq -c | awk '($1==2)' | awk '{print $2}' > SNP_list_for_"$peaks"_"$snps"

####### THIS VERSION IS K27ac ONLY

mkdir AVERAGES
mkdir STAT_"$peaks"_"$snps"_"$dis"_ksonesided


for cells in $cell_types
do
wc -l BINARY_PEAKS_"$peaks"_"$snps"
wc -l HMW_"$peaks"_"$dis"/"$cells"_OUTPUT/"$type"_40_"$dis"/z_scores_final.txt

/lab/corradin_data/merge.plx -a BINARY_PEAKS_"$peaks"_"$snps":4 HMW_"$peaks"_"$dis"/"$cells"_OUTPUT/"$type"_40_"$dis"/z_scores_final.txt:1 | cut -f 5- | sort -u > AVERAGES/zscore_"$peaks"_"$snps"_"$cells"_"$dis"
#wc -l AVERAGES/zscore_"$peaks"_"$snps"_"$cells"_"$dis"

m=$(cat AVERAGES/zscore_"$peaks"_"$snps"_"$cells"_"$dis" | wc -l)
n=$(cat BINARY_PEAKS_"$peaks"_"$snps" | wc -l)

if [ $m == $n ] 
then 
echo ""$cells" counts match"
else 
echo "merge problem for average calculation: "$cells" " 
exit 1
fi

### the numbers above should be equal, exit if merge fail

echo "header	AAAA	peak_bin	R1	R2	R3	R4	R5	R6	R7	R8	R9	R10	R11	R12	R13	R14	R15	R16	R17	R18	R19	R20	R21	R22	R23	R24	R25	R26	R27	R28	R29	R30	R31	R32	R33	R34	R35	R36	R37	R38	R39	R40" > AVERAGES/Averages_"$peaks"_"$snps"_"$cells"_"$dis"

while read s
do
cat AVERAGES/zscore_"$peaks"_"$snps"_"$cells"_"$dis" | grep $s | awk '($2==1)' | awk -v k=$b '{SUM4+=$4} {SUM5+=$5} {SUM6+=$6} {SUM7+=$7} {SUM8+=$8} {SUM9+=$9} {SUM10+=$10} {SUM11+=$11} {SUM12+=$12} {SUM13+=$13} {SUM14+=$14} {SUM15+=$15} {SUM16+=$16} {SUM17+=$17} {SUM18+=$18} {SUM19+=$19} {SUM20+=$20} {SUM21+=$21} {SUM22+=$22} {SUM23+=$23} {SUM24+=$24} {SUM25+=$25} {SUM26+=$26} {SUM27+=$27} {SUM28+=$28} {SUM29+=$29} {SUM30+=$30} {SUM31+=$31} {SUM32+=$32} {SUM33+=$33} {SUM34+=$34} {SUM35+=$35} {SUM36+=$36} {SUM37+=$37} {SUM38+=$38} {SUM39+=$39} {SUM40+=$40} {SUM41+=$41} {SUM42+=$42} {SUM43+=$43} END {print k, $1, "SIG",  SUM4/NR,  SUM5/NR,  SUM6/NR,  SUM7/NR,  SUM8/NR,  SUM9/NR,  SUM10/NR,  SUM11/NR,  SUM12/NR,  SUM13/NR,  SUM14/NR,  SUM15/NR,  SUM16/NR,  SUM17/NR,  SUM18/NR,  SUM19/NR,  SUM20/NR,  SUM21/NR,  SUM22/NR,  SUM23/NR,  SUM24/NR,  SUM25/NR,  SUM26/NR,  SUM27/NR,  SUM28/NR,  SUM29/NR,  SUM30/NR,  SUM31/NR,  SUM32/NR,  SUM33/NR,  SUM34/NR,  SUM35/NR,  SUM36/NR,  SUM37/NR,  SUM38/NR,  SUM39/NR,  SUM40/NR,  SUM41/NR, SUM42/NR, SUM43/NR }' | sed 's/ /\t/g' | sort -dk2 >> AVERAGES/Averages_"$peaks"_"$snps"_"$cells"_"$dis"
cat AVERAGES/zscore_"$peaks"_"$snps"_"$cells"_"$dis" | grep $s | awk '($2==0)' | awk -v k=$b '{SUM4+=$4} {SUM5+=$5} {SUM6+=$6} {SUM7+=$7} {SUM8+=$8} {SUM9+=$9} {SUM10+=$10} {SUM11+=$11} {SUM12+=$12} {SUM13+=$13} {SUM14+=$14} {SUM15+=$15} {SUM16+=$16} {SUM17+=$17} {SUM18+=$18} {SUM19+=$19} {SUM20+=$20} {SUM21+=$21} {SUM22+=$22} {SUM23+=$23} {SUM24+=$24} {SUM25+=$25} {SUM26+=$26} {SUM27+=$27} {SUM28+=$28} {SUM29+=$29} {SUM30+=$30} {SUM31+=$31} {SUM32+=$32} {SUM33+=$33} {SUM34+=$34} {SUM35+=$35} {SUM36+=$36} {SUM37+=$37} {SUM38+=$38} {SUM39+=$39} {SUM40+=$40} {SUM41+=$41} {SUM42+=$42} {SUM43+=$43} END {print k, $1, "NSIG",  SUM4/NR,  SUM5/NR,  SUM6/NR,  SUM7/NR,  SUM8/NR,  SUM9/NR,  SUM10/NR,  SUM11/NR,  SUM12/NR,  SUM13/NR,  SUM14/NR,  SUM15/NR,  SUM16/NR,  SUM17/NR,  SUM18/NR,  SUM19/NR,  SUM20/NR,  SUM21/NR,  SUM22/NR,  SUM23/NR,  SUM24/NR,  SUM25/NR,  SUM26/NR,  SUM27/NR,  SUM28/NR,  SUM29/NR,  SUM30/NR,  SUM31/NR,  SUM32/NR,  SUM33/NR,  SUM34/NR,  SUM35/NR,  SUM36/NR,  SUM37/NR,  SUM38/NR,  SUM39/NR,  SUM40/NR,  SUM41/NR, SUM42/NR, SUM43/NR }' | sed 's/ /\t/g' | sort -dk2 >> AVERAGES/Averages_"$peaks"_"$snps"_"$cells"_"$dis"
done < SNP_list_for_"$peaks"_"$snps"

done


cd AVERAGES 

for cells in $cell_types
do
Rscript "$dir"/Heatmap_Script_Files_PY/Aggregate_ksonesided_paired_greater_tests.r Averages_"$peaks"_"$snps"_"$cells"_"$dis" STAT_"$peaks"_"$snps"_"$cells"_"$dis"
done

#### combine output

echo "HEADER" > LIST_STAT_"$peaks"_"$snps"_"$dis"_ksonesided_temp
for cells in $cell_types
do
echo "STAT_"$peaks"_"$snps"_"$cells"_"$dis"_aggregate_ksonesided_paired_greater_out" >> LIST_STAT_"$peaks"_"$snps"_"$dis"_ksonesided_temp
cat LIST_STAT_"$peaks"_"$snps"_"$dis"_ksonesided_temp | grep -v HEADER > LIST_STAT_"$peaks"_"$snps"_"$dis"_ksonesided

done
cat "$cur_dir"/"$cell_list" | sed 's/ /\t/g' | sed 's/_H3K27ac_hg38//g' | awk '{OFS="\t"}{print "SNP",$0}' > header_"$peaks"_"$snps"_"$dis"_ksonesided
paste $(printf "%s " $(cat LIST_STAT_"$peaks"_"$snps"_"$dis"_ksonesided)) | awk '{OFS="\t"}{print $1,$9,$20,$31,$42,$53,$64,$75,$86,$97,$108,$119,$130,$141,$152,$163,$174,$185,$196,$207,$218,$229,$240,$251,$262,$273,$284,$295,$306,$317,$328,$339,$350,$361,$372,$383,$394}' | cat header_"$peaks"_"$snps"_"$dis"_ksonesided - > ../MATRIX_out_"$peaks"_"$snps"_"$dis"_ksonesided


# else

# echo " ONLY R commands were run"
# cd AVERAGES 

# for cells in $cell_types
# do
# Rscript /lab/corradin_data/permanent/Aggregate_"$r"_paired_greater_tests.r Averages_"$peaks"_"$snps"_"$cells"_"$dis" STAT_"$peaks"_"$snps"_"$cells"_"$dis"
# done

# #### combine output

# echo "HEADER" > LIST_STAT_"$peaks"_"$snps"_"$dis"_"$r"_temp
# for cells in $cell_types
# do
# echo "STAT_"$peaks"_"$snps"_"$cells"_"$dis"_aggregate_"$r"_paired_greater_out" >> LIST_STAT_"$peaks"_"$snps"_"$dis"_"$r"_temp
# cat LIST_STAT_"$peaks"_"$snps"_"$dis"_"$r"_temp | grep -v HEADER > LIST_STAT_"$peaks"_"$snps"_"$dis"_"$r"

# done

# cat ../$cell_list | sed 's/ /\t/g' | sed 's/_H3K27ac_hg38//g' | awk '{OFS="\t"}{print "SNP",$0}' > header_"$peaks"_"$snps"_"$dis"_"$r"
# paste $(printf "%s " $(cat LIST_STAT_"$peaks"_"$snps"_"$dis"_"$r")) | awk '{OFS="\t"}{print $1,$9,$20,$31,$42,$53,$64,$75,$86,$97,$108,$119,$130,$141,$152,$163,$174,$185,$196,$207,$218,$229,$240,$251,$262,$273,$284,$295,$306,$317,$328,$339,$350,$361}' | cat header_"$peaks"_"$snps"_"$dis"_"$r" - > ../MATRIX_out_"$peaks"_"$snps"_"$dis"_"$r"


# fi



