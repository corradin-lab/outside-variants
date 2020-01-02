#!/usr/bin/env python3
#Author: Yanwei Song, Corradin Lab, Whitehead Institute

from datetime import datetime
import os
import sys
import re
import glob
import shutil, subprocess
import processManager_module

working_dir = os.getcwd() #directory from which the script will be executed
SCRIPTS_PATH = working_dir+"/Heatmap_Script_Files_PY/"
params_file = all_chroms = ""


#wigApproach: working with each chromosome wig files
#input: path_to_wig_files, cur_wig_subfolder, chromosomes, num_windows, dist_each_dir, working_dir, SCRIPTS_PATH,num_cores
#output: return sample_chromosomes (chromosomes can find in wig files; list), out_sample_dir
#eg.: wigApproach(path_to_wig_files, cur_wig_subfolder, chromosomes, num_windows, dist_each_dir, working_dir, SCRIPTS_PATH,num_cores)
def wigApproach(path_to_wig_files, cur_wig_subfolder, chromosomes, num_windows, dist_each_dir, working_dir):
    # change working directory to wig subfolder 
    os.chdir(path_to_wig_files+cur_wig_subfolder+"/") 
    # get all .wig file names 
    wig_files = glob.glob("*.wig")
    sample_chromosomes = []
    out_sample_dir = ""
    # working with each chromosomes (we usually have chromosome 1-22)
    for chr in chromosomes:
        valid = 0
        #processing each wig file 
        for wig_file in wig_files:
            #checking if it is wig file (file name usually contains chrxx.wig or chrxx_)
            pattern = f"chr{chr}.wig|chr{chr}_"
            if re.search(pattern, wig_file):
                valid = 1
                #checking if read permissions are set
                if os.access(wig_file, os.R_OK):
                    #if wig file exist and readable, add chr number in
                    sample_chromosomes.append(chr)
                else:
                    print("NOTE:",wig_file,"does not have read permissions. Chromosome", chr, "will be skipped!")
        # if no wig file avaliable for this chromosomes 
        if valid == 0:
            print("WARNING: No wig file for chr"+chr,"found. The chromosome will be skipped!") 

    # if there are wig file avaliable to use 
    if len(sample_chromosomes) > 0: 
        os.chdir(working_dir)
        out_sample_dir = working_dir+"/"+cur_wig_subfolder+"_OUTPUT/wig_"+str(num_windows)+"_"+str(dist_each_dir)+"/"

    #check if output sample directory exist, and overwritten it
    if os.path.exists(out_sample_dir):
        shutil.rmtree(out_sample_dir)
        print(out_sample_dir,"will be overwritten!\n")
    os.makedirs(out_sample_dir)

    #if user want to calculate peak center
    if calc_peak_center == "yes":
        print("Calculating peak medians and splitting",path_to_peak_file,"into individual chromosomes...")
        #calculate medians first, and write in a file named main_peak_file.txt
        calcMedians(out_sample_dir, path_to_peak_file, calc_peak_center)
        #separate main_peak_file.txt by chromosomes
        pickOutChromosomes(calc_peak_center, sample_chromosomes, path_to_peak_file, out_sample_dir)
    else:
        print("Splitting",path_to_peak_file,"into individual chromosomes...")
        #pick out chromosomes from peak file
        pickOutChromosomes(calc_peak_center, sample_chromosomes, path_to_peak_file, out_sample_dir)
    print("Done\n")

    print ("Scheduling parallel processing jobs and running HeatmapWindowMedians on",cur_wig_subfolder,"...\n\n")

    return sample_chromosomes, out_sample_dir



#bigwigApproach: working with bigwig files, if there are bigwig files in this subfolder, then return out_sample_dir. 
#input: working_dir, cur_wig_subfolder, num_windows, dist_each_dir
#output: return out_sample_dir, bigwig
#eg.: bigwigApproach(working_dir, cur_wig_subfolder, num_windows, dist_each_dir, chromosomes, path_to_wig_files)
############right now, use convert bigwig to wig file approach#######
def bigwigApproach(working_dir, cur_wig_subfolder, num_windows, dist_each_dir, chromosomes, path_to_wig_files):
    sample_bigwig = {}
    # change working directory to bigwig subfolder 
    os.chdir(path_to_wig_files) 
    # get all .bw file names 
    bw_files = glob.glob(cur_wig_subfolder+".*")
    bigwig = ""
    if len(bw_files) > 0:
        out_sample_dir = working_dir+"/"+cur_wig_subfolder+"_OUTPUT/bigwig_"+str(num_windows)+"_"+str(dist_each_dir)+"/"
        for bw_file in bw_files:
            pattern = f"bw|bigWig|bigwig"
            if re.search(pattern, bw_file) and os.path.getsize(bw_file) > 0:
                if os.access(bw_file, os.R_OK):
                    sample_bigwig[bw_file] = os.path.getsize(bw_file)
                else:
                    print("NOTE:",bw_file,"does not have read permissions. Please check file permissions!")
    else: 
        sys.exit("There's no bigwig file found in "+path_to_wig_files+". Please check it again!!!")
        
    #check if output sample directory exist, and overwritten it 
    if os.path.exists(out_sample_dir):
        shutil.rmtree(out_sample_dir)
        print(out_sample_dir,"will be overwritten!\n")
    os.makedirs(out_sample_dir)

    #check if bigwig file read permissions are set
    if len(sample_bigwig) > 1:
        bigwig = list(sample_bigwig.keys())[0] #first bigwig file will be used 
        print("WARNING: there are multiple bigwig files for",cur_wig_subfolder,"only", bigwig, "will be used in this analysis.")
        #convert bigwig to wig base on chromosomes user provided
    elif len(sample_bigwig) == 1:
        bigwig = list(sample_bigwig.keys())[0]
    else:
        sys.exit("There's no bigwig file avaliable to use",path_to_wig_files,". Please check it again!!!")

    #if user want to calculate peak center 
    if calc_peak_center == "yes":
        print("Calculating peak medians",path_to_peak_file+"...")
        calcMedians(out_sample_dir, path_to_peak_file, calc_peak_center) #only need to calculate median for each peaks in bigwig file 
        #separate main_peak_file.txt by chromosomes
        pickOutChromosomes(calc_peak_center, chromosomes, path_to_peak_file, out_sample_dir)
    else: 
        print("Splitting",path_to_peak_file,"into individual chromosomes...")
        #pick out chromosomes from peak file
        pickOutChromosomes(calc_peak_center, chromosomes, path_to_peak_file, out_sample_dir)
        # print("No peak medians need to calculate, use peak file",path_to_peak_file,"directly.")
    print("Done\n")

    print ("Scheduling parallel processing jobs and running HeatmapWindowMedians on",cur_wig_subfolder,"...\n\n")

    return out_sample_dir, bigwig






#calcMedians: if calculate peak center is yes, calculate center (median) for each peak, write into main_peak_file.txt
#input: path_to_peak_file, out_sample_dir, sample_chromosomes
#output: write results in OUTPUT files 
#eg.: calcMedians(out_sample_dir, path_to_peak_file, calc_peak_center)
def calcMedians(out_sample_dir, path_to_peak_file, calc_peak_center):
    try:
        OUTPUT = open(out_sample_dir+"main_peak_file.txt", "w")
    except IOError:
        sys.exit("Problem calculating peak centers\n")
    try:
        INPUT = open(path_to_peak_file)
    except IOError:
        sys.exit("Problem calculating peak centers\n")
    else: 
        with INPUT:
            for line in INPUT:
                line = line.strip()
                if re.match('chr[^\t\s]+[\t\s]+\d+[\t\s]+\d+.*', line):
                    if calc_peak_center == "yes":
                        l = line.split()
                        if (int(l[1]) + int(l[2])) % 2 == 0:
                            peak_center = int((int(l[1]) + int(l[2]))/2) 
                            OUTPUT.write(l[0]+"\t"+str(peak_center)+"\n")
                        else:
                            peak_center = eval("(int(l[1]) + int(l[2])) / 2")
                            OUTPUT.write(l[0]+"\t"+str(peak_center)+"\n")
                else: 
                    sys.exit("The peak file is improperly formatted.The peak file must be in this format: chromosome< TAB or space>< midpoint coord>... . The program will exit.\n")
    OUTPUT.close()



#pickOutChromosomes: separate peaks file or main_peak_file.txt file based on chromosomes. 
#input: path_to_peak_file, path_to_peak_file, sample_chromosomes
#output: nothing return, separate peak file based on chr and write into out sample folder 
#eg.: pickOutChromosomes(calc_peak_center, sample_chromosomes, path_to_peak_file, out_sample_dir)
def pickOutChromosomes(calc_peak_center, sample_chromosomes, path_to_peak_file, out_sample_dir):
    peakFile = ""
    if calc_peak_center == "yes":
        peakFile = out_sample_dir+"main_peak_file.txt"
    else:
        peakFile = path_to_peak_file
    for chr in sample_chromosomes:
        os.system("awk '$1==\"chr"+chr+"\"' "+peakFile+" > "+out_sample_dir+"/chr"+chr+".txt") #or die "Can't extract chr$chr peaks from main_peak_file"
        os.system("sort "+out_sample_dir+"chr"+chr+".txt -n -k 2 -u > "+out_sample_dir+"tmp"+chr+".txt")
        os.system("mv "+out_sample_dir+"tmp"+chr+".txt "+out_sample_dir+"chr"+chr+".txt")
    os.system("rm "+out_sample_dir+"main_peak_file.txt")






#getParams: read commandline parameters and take values from param file
#input: nothing, will read input arguments in the function 
#output: return chromosomes, path_to_peak_file, path_to_wig_files, wig_subfolders, calc_peak_center, num_windows, dist_each_dir, fileType
#eg.: getParams()
def getParams():
    #if user do not provide param file name and proper number of core, stop and give information 
    if len(sys.argv) < 3:
        sys.exit("Usage: ./getPeakWindowMedians params_file.txt num_cores\n")
    params_file = sys.argv[1]
    num_cores = int(sys.argv[2])
    if num_cores < 0 or num_cores>40:
        sys.exit("You must specify between 1 and 40 cores\n")

    chromosomes = []
    wig_subfolders = []
    path_to_peak_file = ""
    fileType = ""
    calc_peak_center = ""
    path_to_wig_files = ""
    num_windows = 0
    dist_each_dir = 0
    ##read params file parameters, if file can't read then stop and provide error information 
    try:
        PARAMS_FILE = open(params_file, 'r+')
    except IOError:
        sys.exit("Invalid parameters file\n")
    else:
        #open param file and extract: chromosomes, path to peak file, path to wig files, wig subfolders, calculate peak center, number of windows and distant from peak in each directory
        with PARAMS_FILE:
            for l in PARAMS_FILE:
                line = l.strip().split(":")
                if line[0] == "chromosomes":
                    if "," in line[1]:
                        chromosomes = line[1].split(",")
                    else:
                        chromosomes.append(line[1])
                elif line[0] == "absolute_path_to_peak_file":
                    try: 
                        peakF = open(line[1])
                        peakF.close()
                    except FileNotFoundError:
                        sys.exit("peak file",line[1], "doesn't exist")
                    else: 
                        path_to_peak_file = line[1]
                elif line[0] == "absolute_path_to_sample_dirs":
                    if line[1].endswith("/"):
                        path_to_wig_files = line[1]
                    else:
                        path_to_wig_files = line[1]+"/"
                    # if not os.path.exists(path_to_wig_files):
                    #     sys.exit("wig directory", line[1], "doesn't exist")
                elif line[0] == "sample_dirs":
                    if "," in line[1]:
                        wig_subfolders = line[1].split(",")
                        new_wig_subfolders = []
                        for i in wig_subfolders:
                            if os.path.exists(path_to_wig_files+i+"/"):
                                new_wig_subfolders.append(i)
                            else:
                                print("WARNING: sample", i, "doesn't exist and will be ignored")
                        wig_subfolders = new_wig_subfolders
                    else:
                        if os.path.exists(path_to_wig_files+line[1]+"/"):
                            wig_subfolders.append(line[1])
                        else:
                            print(path_to_wig_files+line[1])
                            sys.exit("sample directory", path_to_wig_files+line[1]+'/', "doesn't exist")
                elif line[0] == "calculate_peak_center?":
                    calc_peak_center = line[1]
                elif line[0] == "num_windows":
                    num_windows = int(line[1])
                elif line[0] == "dist_from_peak_in_each_direction":
                    dist_each_dir = int(line[1])
                elif line[0] == "file_type":
                    fileType = line[1].lower()
                    if fileType == "wig":
                        if not os.path.exists(path_to_wig_files):
                            sys.exit("wig directory", line[1], "doesn't exist")
                else:
                    print("Not recognized as a parameter line:", l, "\n")

    if len(chromosomes) <= 0:
        sys.exit("You must specify at least 1 chromosome\n")
    elif path_to_peak_file == "":
        sys.exit("You must provide the full path to a valid regions file\n")
    elif path_to_wig_files == "":
        sys.exit("You must provide the full path to sample directories (including all forward slashes)\n")
    elif len(wig_subfolders) <= 0:
        sys.exit("You must specify at least one sample directory (no forward slashes)\n")
    elif calc_peak_center == "":
        sys.exit("You must specify if peak medians will need to be calculated\n")
    elif num_windows <= 0:
        sys.exit("You must provide number of peak windows\n")
    elif dist_each_dir <= 0:
        sys.exit("You must specify the distance in each direction of the peak the median signal will need to be calculated for\n")
    elif fileType == "":
        sys.exit("You must specify the which type of sample file you want to use (wig or bigwig file)\n")
    print("Finished params validation\n")
    return chromosomes, path_to_peak_file, path_to_wig_files, wig_subfolders, calc_peak_center, num_windows, dist_each_dir, fileType, num_cores
            



#z_score: calculate z-score for reformated out sample file and write results in z_scores_filnal.txt for each sample in each wig subfolder
#input: whole path to the FINAL_ALL_CHR.txt file (usually in out sample file folder), name of out sample folder, path to wig subfolders
#output: no return, write results in an output file 
#eg.: z_score(out_sample_dir+"FINAL_ALL_CHR.txt", cur_wig_subfolder, out_sample_dir, SCRIPTS_PATH, num_windows)
def z_score(file_to_z_score, sample, output_dir, SCRIPTS_PATH, num_windows): 
    OUTPUT = open(output_dir+'z_scores_final.txt', 'w')
    max = num_windows + 1
    # call Rscript to generate an temp.txt file contains z-scores
    command = "R --vanilla "+file_to_z_score+" < "+SCRIPTS_PATH+"z_score.R"
    os.system(command)
    OUTPUT.write("chr_coord\t")

    #give header line for in output file (each window sizes represent in differnet column), and add z-score file in output file
    for i in range(1,num_windows,1):
        OUTPUT.write(sample+"_"+str(i)+"\t")
    OUTPUT.write(sample+'_'+str(num_windows)+'\n')
    OUTPUT.close()
    command2 = "cut -f1-"+str(max)+" "+output_dir+"temp.txt >> "+output_dir+"z_scores_final.txt"
    os.system(command2)






chromosomes, path_to_peak_file, path_to_wig_files, wig_subfolders, calc_peak_center, num_windows, dist_each_dir, fileType, num_cores = getParams()

print("-----------Main program started on: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+"-------------\n")


# working with each wig subfolder 
for cur_wig_subfolder in wig_subfolders: 
    if fileType == "wig": # for wig files
        sample_chromosomes,out_sample_dir = wigApproach(path_to_wig_files, cur_wig_subfolder, chromosomes, num_windows, dist_each_dir, working_dir)
        bw_file = None
        # if the number of out sample chr less than user want to run
        if len(sample_chromosomes) <= len(chromosomes) and len(sample_chromosomes) > 0:
            sample_chromosomes = sample_chromosomes
        processManager_module.processManager(sample_chromosomes, cur_wig_subfolder, path_to_wig_files, out_sample_dir, num_windows, dist_each_dir, num_cores, bw_file, fileType)

    elif fileType == "bigwig": #for bigwig files
        out_sample_dir, bw_file = bigwigApproach(working_dir, cur_wig_subfolder, num_windows, dist_each_dir, chromosomes, path_to_wig_files)

        processManager_module.processManager(chromosomes, cur_wig_subfolder, path_to_wig_files, out_sample_dir, num_windows, dist_each_dir, num_cores, bw_file, fileType, SCRIPTS_PATH)
        os.system('rm '+out_sample_dir+'chr*.wig') #remove wig files (only works for bigwig approach)

    # using shell command to generate duplicate removed final file
    os.system('cat '+out_sample_dir+'FINAL*_NO_DUPS > '+out_sample_dir+'FINAL_ALL_CHR.txt')
    os.system('rm '+out_sample_dir+'FINAL*_NO_DUPS')
    #######Calculate z-score########
    print("Z-scoring...")
    os.chdir(out_sample_dir)
    z_score(out_sample_dir+"FINAL_ALL_CHR.txt", cur_wig_subfolder, out_sample_dir, SCRIPTS_PATH, num_windows)
    print("Done\n")
    os.chdir(working_dir)
    ##Clean up###


print("-----------Main program finished on: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+"-------------\n")







