#!/usr/bin/env python3
#Author: Yanwei Song, Corradin Lab, Whitehead Institute
#eg.: getWindowMedians.py 

import re
import sys
import os
import subprocess
import glob
from decimal import Decimal


#main script 
def getWindowMedians(args):
    fileType, full_path_to_sample_wigs, output_dir, num_windows, dist_each_dir, chr, bw_file, SCRIPTS_PATH = args
    pnum = os.getpid()
    wig_signal_vals = {}
    peak_windows = {}

    print("Process", str(pnum)+": Starting chromosome", chr+"...")

    if fileType == "bigwig": #convert bw file to wig file 
        step = 10 #right now the default step is 10, we can ask user provide the step they want to use 
        convBWtoWig(bw_file, chr, step, output_dir, SCRIPTS_PATH)
        full_path_to_sample_wigs = output_dir

    #check if wig file and peak file can open successfully. 
    wig_file = validateChromWig(full_path_to_sample_wigs, chr)
    ##check wig file and get wig file info
    peakfile_open, CHR_PEAK_FILE, wig_open, CHR_WIG_FILE = wigFileChecking(chr,output_dir,wig_file) 
    #checks to make sure wig file is properly formatted
    success, type_wig, wig_step, wig_start_coord = getWigFileInfo(wig_file, chr)

    if wig_open == peakfile_open == success == "true": #proceed only of both the wig and peak files are open 
        try:
            OUTFILE = open(output_dir+"chr"+chr+"_OUT", "w")
        except IOError:
            print("Can't open peaks file:",output_dir+"chr"+chr+".txt\n")
        #store peaks 
        peaks = populatePeakList(CHR_PEAK_FILE)
        window_size = dist_each_dir*2/num_windows
        #read wig file line by line, compare with peak file, calculate medians in each window of each peak and write it in output file 
        wigFileReading(wig_file, CHR_WIG_FILE, chr, peaks, dist_each_dir, window_size, num_windows, OUTFILE, wig_signal_vals, peak_windows)
        #finish looping through wig file 

        for peak in peaks: 
            for i in range(1, num_windows+1):
                OUTFILE.write(str(peak)+"\tNA\n")

        OUTFILE.close()
        wig_signal_vals = {}
        peak_windows = {}
    #finish checking if both wig and peak file are open


    #######Convert to right format & calculate Z-score########
    #process out (_OUT) file 
    cur_file = "chr"+chr+"_OUT"
    print("Converting",cur_file,"to the right format...")
    removeDuplicates(cur_file, output_dir)
    convertToRightFormat(cur_file, output_dir, num_windows)
    print("Process", str(pnum)+": Done with chromosome", chr+"...\n")

    #finish checking if both wig and peak file are open 
#finish looping through chromosomes 



#convBWtoWig: convert bigwig file to wig files. 
#input: chr number, and step, bigwig file name, path to output directory 
#output: no return, write wig files out
#eg.: convBWtoWig(bw_file, chr, step, out_sample_dir)
def convBWtoWig(bw_file, chr, step, out_sample_dir, SCRIPTS_PATH):
    wig = out_sample_dir+"chr"+chr+".wig"
    tempF = out_sample_dir+"temp"+chr+".bedgraph"
    os.system("bigWigToWig "+bw_file+" "+tempF+" -chrom=chr"+chr)
    os.system("perl "+SCRIPTS_PATH+"bedgraph_to_wig.pl --bedgraph "+tempF+" --wig "+wig+" --step "+str(step))
    os.system("rm "+tempF)


#validateChromWig: check if wig file of chromosome number is valid 
#input: full_path_to_sample_wigs, chr
#output: return a wig file name. 
#e.g.: validateChromWig(full_path_to_sample_wigs, chr)
def validateChromWig(full_path_to_sample_wigs, chr):
    pattern = f"*chr{chr}[.(_.)]*wig"
    wig_files = glob.glob(full_path_to_sample_wigs+pattern) # get all wig files in this directory of this chr
    if len(wig_files) == 0:
        print("Could not find any chr"+chr,"wig files in", full_path_to_sample_wigs, "in this format: *chr",chr+".wig* or *chr"+chr+"_*.wig. Skipping chr"+chr)
        return 0 
    elif len(wig_files) > 1:
        new_wig_files = []
        for i in wig_files:
            if os.stat(i).st_size != 0:
                new_wig_files.append(i)
        if len(new_wig_files) > 1:
            print("WARNING: there are",len(new_wig_files), "chr"+chr,"wig files find in",full_path_to_sample_wigs+",",new_wig_files[0],"will be used in this analysis.")
            return new_wig_files[0]
        elif len(new_wig_files) == 1:
            print("Only one nonempty wig file of chr"+chr,"found,",new_wig_files,"will be used in this analysis.")
            return new_wig_files[0]
        else:
            print("WORNING: all chr"+chr,"wig files are empty, no wig file available to use. Skipping chr"+chr)
            return 0
    else:
        if os.stat(wig_files[0]).st_size == 0:
            print("WORNING:",wig_files[0],"appears to be empty, no wig file available to use. Skipping chr"+chr)
            return 0
        elif os.stat(wig_files[0]).st_size > 0:
            print("Only one wig file of chr"+chr,"found,",wig_files[0],"will be used in this analysis.")
            return wig_files[0]



#wigFileChecking: check wig file format and if can open properly for each chromosomes and if peak can open successfully. 
#input: pnum (number of processors; int), chr (individual chromosome number; int), full_path_to_sample_wigs (/path/to/wigfile/sampelType/; list), output_dir (/path/to/out/sample/; str), wig_file (wig file name, get from validateChromWig; str)
#output: return peakfile_open (true or false; str), wig_open (true or false; str)
#e.g.: peakfile_open, CHR_PEAK_FILE, wig_open, CHR_WIG_FILE = wigFileChecking(pnum, chr, full_path_to_sample_wigs,output_dir,wig_file) 
def wigFileChecking(chr,output_dir,wig_file):
    wig_open = "true"   
    peakfile_open = "true"
    try:
        CHR_WIG_FILE = open(wig_file, "r")
    except IOError:
        print("Can't open",wig_file,"\n")
        wig_open = "false"
    try:
        CHR_PEAK_FILE = open(output_dir+"chr"+chr+".txt", "r")
    except IOError:
        print("Can't open",output_dir+"chr"+chr+".txt\n")
        peakfile_open = "false"
    return peakfile_open, CHR_PEAK_FILE, wig_open, CHR_WIG_FILE



#getWigFileInfo: checking if wig file can open properly, identify wig file type based on header line, if there are any variable missed based different wig file type.
#input: wig_file, chr
#output: return success(true/false), wig_type (variableStep/fixedStep), step (an integer), start (an integer). 
#e.g.: getWigFileInfo(wig_file, chr)
def getWigFileInfo(wig_file, chr):
    wig_type = ""
    step = ""
    start = ""
    chrom = ""
    success = ""

    wig_header_line = subprocess.check_output("grep -E 'fixedStep|variableStep' "+wig_file, shell = True).decode("utf-8").rstrip()
    wig_header_line_charas = wig_header_line.split()
    #if chromosome number include in header 
    if re.match("chrom=(\S+)", wig_header_line_charas[1]):
        if wig_header_line_charas[1].split("=")[1] == "chr"+chr:
            chrom = wig_header_line_charas[1].split("=")[1]
    else:
        success = "false"
        print(wig_file, "does not contain chr"+chr, "information. Chromosome", chr, "will be skipped\n")
    
    if re.match("fixedStep", wig_header_line):
        wig_type = "fixedStep"
        if re.match("start=(\d+)", wig_header_line_charas[2]):
            start = wig_header_line_charas[2].split("=")[1]
        if re.match("step=(\d+)", wig_header_line_charas[3]):
            step = wig_header_line_charas[3].split("=")[1]
    elif re.match("variableStep", wig_header_line):
        wig_type = "variableStep"
        step = "NA"
        start = "NA"
    else:
        print(wig_file,"is malformed. Chromosome",chr,"will be skipped\n")
        success = "false"
    
    if wig_type != "" and step != "" and start != "" and chrom != "":
        success = "true"
    else:
        print(wig_file,"is malformed. Chromosome",chr,"will be skipped\n")
        success = "false"
    return success, wig_type, step, start



#populatePeakList: extract start position of each peak in chromosome peaks file, and check duplication of them.
#input: CHR_PEAK_FILE (file handle of peak file for each chromosome)
#output: return a list, contains no_dup peaks from CHR_PEAK_FILE. 
#e.g.: populatePeakList(CHR_PEAK_FILE)
def populatePeakList(CHR_PEAK_FILE):
    peaks = []
    existing_peaks = {}
    with CHR_PEAK_FILE:
        for peakline in CHR_PEAK_FILE:
            columns = peakline.strip().split()
            if columns[1] in existing_peaks.keys():
                existing_peaks[columns[1]] += 1
            else: 
                try:
                    peaks.append(int(columns[1]))
                except ValueError:
                    peaks.append(float(columns[1]))
                existing_peaks[columns[1]] = 1
    return peaks



#wigFileReading: use getWigFileInfo to get wig file information and read it line by line, compare wig file with peak file and use calcMediansInWindow to calculate median for each window and write in an output file. 
# check wig file format and if can open properly for each chromosomes and if peak can open successfully. 
#input: wig_file (wig file name, get from validateChromWig; str), CHR_WIG_FILE (filehandle of wig file), chr (individual chromosome number; int), peaks (dup removed peaks; list), dist_each_dir (distance from each direction; int), window_size (get from distance of each direction and number of windowns, dist_each_dir*2/num_windows; int) 
# pnum (number of processors; int), chr (individual chromosome number; int), full_path_to_sample_wigs (/path/to/wigfile/sampelType/; list), output_dir (/path/to/out/sample/; str)
#output: return peakfile_open (true or false; str), wig_open (true or false; str)
#e.g.: wigFileReading(wig_file, CHR_WIG_FILE, chr, peaks, dist_each_dir, window_size, num_windows, OUTFILE, wig_signal_vals, peak_windows)
def wigFileReading(wig_file, CHR_WIG_FILE, chr, peaks, dist_each_dir, window_size, num_windows, OUTFILE, wig_signal_vals, peak_windows):
    num_wig_steps = 0
    flag = 0
    peak = 0
    wig_coord = 0
    count = 0

    size = subprocess.check_output("wc "+wig_file, shell = True).decode("utf-8").rstrip().split()[0]

    success, type_wig, wig_step, wig_start_coord = getWigFileInfo(wig_file, chr) #checks to make sure wig file is properly formatted
    wig_start_coord = int(wig_start_coord)
    wig_step = int(wig_step)
    with CHR_WIG_FILE:
        for wig_line in CHR_WIG_FILE:
            count += 1 
            wig_line = wig_line.strip()
            #if data line, then proceed (must start with integer)
            if wig_line[0].isdigit():
                if type_wig == "fixedStep":
                    wig_coord = wig_start_coord+num_wig_steps*wig_step-1
                    # wig_signal_val = float(wig_line)
                    try:
                        wig_signal_val = (int(wig_line))
                    except ValueError:
                        wig_signal_val = (float(wig_line))
                else: #variable step 
                    wig_liner = wig_line.split()
                    wig_coord = int(wig_liner[0])-1
                    wig_signal_val = float(wig_liner[1])
                    if re.search("E", wig_signal_val): #convert scientific notation to integer
                        wig_signal_val = float(wig_signal_val)
                num_p_window = 1
                #remove peaks that occur before current coordinate analyzed 
                if flag > 0:
                    del peaks[:flag]
                    flag = 0

                for peak in peaks:
                    if not peak_windows.get(peak):
                        peak_windows[peak] = 1

                    num_p_window = peak_windows[peak]
                    peak_start = peak - dist_each_dir + ((num_p_window - 1) * window_size)
                    peak_end = peak - dist_each_dir + (num_p_window * window_size)

                    if count == int(size):
                        if wig_coord < peak_end and wig_coord >= peak_start:
                            if peak in wig_signal_vals.keys():
                                wig_signal_vals[peak].append(wig_signal_val)
                            else:
                                wig_signal_vals[peak] = [wig_signal_val]
                        if peak in wig_signal_vals.keys():
                            calcMediansInWindow(wig_signal_vals[peak], peak, peak_start, peak_end, OUTFILE)
                            if peak in wig_signal_vals.keys() and len(wig_signal_vals[peak]) > 0:
                                peaks.pop(0)

                        wig_signal_vals[peak] = None
                        peak_end2 = peak_end
                        peak_start2 = peak_start

                        while peak_windows[peak] < num_windows:
                            calcMediansInWindow("NA", peak, peak_start2, peak_end2, OUTFILE)
                            peak_windows[peak] += 1
                            peak_start2 = peak - dist_each_dir + (peak_windows[peak] - 1) * window_size
                            peak_end2 = peak - dist_each_dir + (peak_windows[peak] * window_size)
                        if peak_windows[peak] >= num_windows:
                            if wig_coord >= peak_end2 and peak_windows[peak] == num_windows: 
                                final_start = peak - dist_each_dir + ((peak_windows[peak] - 1) * window_size)
                                final_end = peak - dist_each_dir + (peak_windows[peak] * window_size)
                                calcMediansInWindow("NA", peak, final_start, final_end, OUTFILE) 
                            flag += 1 
                            del wig_signal_vals[peak]
                            del peak_windows[peak]
                    #end of eof

                    # if wig coord falls within given peak window, store wig signal 
                    elif wig_coord < peak_end and wig_coord >= peak_start: 
                        if peak in wig_signal_vals.keys():
                            wig_signal_vals[peak].append(wig_signal_val)
                        else:
                            wig_signal_vals[peak] = [wig_signal_val]
                    # coordinate exceeds allowed range 
                    elif wig_coord >= peak_end:
                        if peak in wig_signal_vals.keys():
                            calcMediansInWindow(wig_signal_vals[peak], peak, peak_start, peak_end, OUTFILE)
                        num_p_window = peak_windows[peak]
                        wig_signal_vals[peak] = None
                        peak_end2 = peak_end
                        peak_start2 = peak_start

                        while wig_coord >= peak_end2 and peak_windows[peak] < num_windows:
                            calcMediansInWindow("NA", peak, peak_start2, peak_end2, OUTFILE)
                            peak_windows[peak] += 1 
                            peak_start2 = peak - dist_each_dir + ((peak_windows[peak] - 1) * window_size)
                            peak_end2 = peak - dist_each_dir + (peak_windows[peak] * window_size)

                        if peak_windows[peak] <= num_windows and wig_coord < peak_end2:
                            if wig_signal_vals.get(peak):
                                wig_signal_vals[peak].append(wig_signal_val)
                            else:
                                wig_signal_vals[peak] = [wig_signal_val]
                        else: #begin else  
                            if wig_coord >= peak_end2 and peak_windows[peak] == num_windows:
                                final_end = peak - dist_each_dir + (peak_windows[peak] * window_size)
                                final_start = peak - dist_each_dir + ((peak_windows[peak] - 1) * window_size)
                                calcMediansInWindow("NA", peak, final_start, final_end, OUTFILE)
                            flag += 1 
                            del wig_signal_vals[peak]
                            del peak_windows[peak]
                        #end else 

                    #wig coordinate is smaller than any of the coordinates in the current window of the current peak. We can skip all subsequent peaks because that will also be true for those peaks. 
                    else: 
                        break
                #finish looping through peaks 
                num_wig_steps += 1
            else: #if not data line, don't do anything 
                pass
    #finish looping through wig file 



#calcMediansInWindow: calculate medians and write in outfile
#input: wig_signal_vals[peak], peak, peak_start, peak_end
#output: nothing return, write peak, range (start and end) to output files
#e.g.: calcMediansInWindow(wig_signal_vals[peak], peak, peak_start, peak_end, OUTFILE)
def calcMediansInWindow(pole, peak, peak_start, peak_end, OUTFILE):
    if not isinstance(pole, list):
        pole = [pole]
    ret = 0
    if len(pole) >= 0 and pole[0] != "NA":
        pole.sort()
        if len(pole) % 2 == 1:
            ret = pole[int(((len(pole) + 1) / 2)) - 1]
        else:
            ret = (pole[int((len(pole) / 2) - 1)] + pole[int(len(pole) / 2)]) / 2
    else:
        ret = "NA"
    OUTFILE.write(str(peak)+"\t"+str(ret)+"\trange:\t"+str(peak_start)+"\tand\t"+str(peak_end)+"\n")




#removeDuplicates: if 1st, 4th (range start) and 6th (range end) same in the same time, then we defined this line is duplicated and remove duplication. 
#input: cur_file, out_sample_dir
#output: nothing return, output file is duplicate removed file 
#eg.: removeDuplicates(cur_file, out_sample_dir )
def removeDuplicates(cur_file, out_sample_dir):
    repeats = {}
    OUTPUT = open(out_sample_dir+cur_file+"_NO_DUPS", 'w+')
    with open(out_sample_dir+cur_file) as INPUT:
        for line in INPUT: 
            line = line.strip()
            l = line.split()
            if len(l) > 2:
                k = l[0]+" and "+l[3]+" and "+l[5]
            else:
                k = l[0]
            if k not in repeats.keys():
                repeats[k] = 1
                OUTPUT.write("\t".join(l)+"\n")
    OUTPUT.close()
    os.system("rm "+out_sample_dir+cur_file)



#convertToRightFormat: 
#input:
#output: 
#eg.: convertToRightFormat(cur_file, out_sample_dir, num_windows)
def convertToRightFormat(input_file, output_dir, num_windows):
    input_file = input_file+'_NO_DUPS'
    peaks = {}
    sizes = {}
    OUTPUT = open(output_dir+"FINAL_"+input_file, 'w')

    name = input_file.split("_")
    chr_name = name[0]

    with open(output_dir+input_file) as INPUT:
        for line in INPUT:
            line = line.strip()
            columns = line.split()

            if columns[0] not in peaks.keys():
                sizes[columns[0]] = 1
            else:
                sizes[columns[0]] = sizes[columns[0]]+1
            size = sizes[columns[0]]
            if size <= num_windows:
                if columns[0] in peaks.keys():
                    peaks[columns[0]].append(columns[1])
                else:
                    peaks[columns[0]] = [columns[1]]
            if size == num_windows:
                arr = peaks[columns[0]]
                OUTPUT.write(chr_name+"\t"+columns[0]+"\t")
                OUTPUT.write("\t".join(arr)+"\n")
                #new line of code
                del peaks[columns[0]]
    for k in peaks.keys():
        arr = peaks[k]
        OUTPUT.write(chr_name+'\t'+k+'\t')
        if num_windows >= 0:
            OUTPUT.write("\t".join(arr)+"\n")
    OUTPUT.close()









