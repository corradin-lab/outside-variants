#!/usr/bin/env python3
#Author: Yanwei Song, Corradin Lab, Whitehead Institute
#eg.: processManager.py all_chroms cur_wig_subfolder path_to_wig_files out_sample_dir num_windows dist_each_dir num_cores

import os
import sys
import operator
import getWindowMedians_module 
from multiprocessing import Pool
import time




#ChromChecking: checking if peak file contains that chromosomes, and calculate the number of processors will needed. 
#input: all_chroms(chromosome names; list), cur_wig_subfolder (sample type; list), path_to_wig_files(/path/to/wig/file/; str),output_dir (/path/to/out/sample/; str), num_windows (number of windows want to run; int), dist_each_dir (distance to each directions; int), num_cores(number of cores; int)
#output:  
#e.g.: ChromChecking(chroms_to_run, cur_wig_subfolder, path_to_wig_files, output_dir,num_windows, dist_each_dir, num_cores)
# def ChromChecking(chroms_to_run, cur_wig_subfolder, path_to_wig_files, output_dir, num_windows, dist_each_dir, num_cores):
def ChromChecking(chroms_to_run, output_dir, num_cores):
    file_sizes = {}
    chr_count = len(chroms_to_run)
    num_processors = min(num_cores,40,chr_count)
    for chr in chroms_to_run:
        #if peak file exist then get file sizes 
        if os.path.getsize(output_dir+"chr"+chr+".txt") == 0: 
            print ("WARNING: The peak file for chr"+chr," appears to be empty. The chromosome will not be processed\n")
        elif os.path.getsize(output_dir+"chr"+chr+".txt") > 0:
            file_sizes[chr] = os.path.getsize(output_dir+"chr"+chr+".txt")
        else: 
            file_sizes[chr] = 100
    num_valid_chroms = len(file_sizes) ##This will remove all chromosomes not found in peak file
    #if number of valid chromosomes is changed to smallest value here, number of processors will change 
    num_processors = min(num_cores, 40, num_valid_chroms)

    chr_k = list(file_sizes.keys())

    ##Approach 1: call ProcessorSchedule directly in this function
    # return num_processors, sorted_chr_values, sorted_chr_keys, chr_count
    # ProcessorSchedule(num_processors,sorted_chr_values, sorted_chr_keys, chr_count, cur_wig_subfolder, path_to_wig_files, output_dir, num_windows, dist_each_dir)

    ##Approach 2: use multiprocessing function, call MultiProcessManage directly here
    # return num_processors, sorted_chr_keys
    return num_processors, chr_k
    # MultiProcessManage(num_processors, sorted_chr_keys, cur_wig_subfolder, path_to_wig_files, output_dir, num_windows, dist_each_dir)



#MultiProcessManage: schedule processors
#input: num_processors(number of processor may needed; int), sorted_chr_values(sizes of peak files, sorted based on chr number (decending); list), sorted_chr_keys (chromosome number, decending sorted; list), chr_count (number of chromosomes need to run; int), cur_wig_subfolder (sample type; list), path_to_wig_files(/path/to/wig/file/; str), output_dir (/path/to/out/sample/; str), num_windows (number of windows want to run; int), dist_each_dir (distance to each directions; int), fileType (wig file or bigwig file; str)
#output: 
#e.g.: MultiProcessManage(num_processors, sorted_chr_values, sorted_chr_keys,chr_count,cur_wig_subfolder, path_to_wig_files, output_dir, num_windows, dist_each_dir)
def MultiProcessManage(fileType, num_processors, sorted_chr_keys, cur_wig_subfolder, path_to_wig_files, output_dir, num_windows, dist_each_dir, bw_file, SCRIPTS_PATH):
    #create a pool contains the number of processors can use 
    po = Pool(processes = num_processors)
    full_path_to_sample_wigs = path_to_wig_files+cur_wig_subfolder+"/"
    print(full_path_to_sample_wigs)
    #add all works in pool
    for chr in sorted_chr_keys:
        args = [fileType, full_path_to_sample_wigs, output_dir, num_windows, dist_each_dir, chr, bw_file, SCRIPTS_PATH] #array of chroms assigned to pid, pid, cur_wig_subfolder, path_to_wig_files, output directory        
        po.apply_async(runWindowMedians, args=(args,))
    #finish looping through chromosomes

    po.close()
    po.join()
    print("Subprocesses done.\n\n")



#runWindowMedians: use function in getWindowMedians module 
#input: args
#output: nothing return 
#e.g.: runWindowMedians(args)
def runWindowMedians(args):
    getWindowMedians_module.getWindowMedians(args) # use functions in getWindowMedians_module.pl 




# def processManager(chroms_to_run, cur_wig_subfolder, path_to_wig_files, output_dir, num_windows, dist_each_dir, SCRIPTS_PATH, num_cores, bw_file, fileType):
def processManager(chroms_to_run, cur_wig_subfolder, path_to_wig_files, output_dir, num_windows, dist_each_dir, num_cores, bw_file, fileType, SCRIPTS_PATH):
    #check if peak file contains this chromosome
    num_processors, sorted_chr_keys = ChromChecking(chroms_to_run, output_dir, num_cores)

    #run multi-jobs parallely 
    MultiProcessManage(fileType, num_processors, sorted_chr_keys, cur_wig_subfolder, path_to_wig_files, output_dir, num_windows, dist_each_dir, bw_file, SCRIPTS_PATH)







#ProcessorSchedule: schedule processors
#input: num_processors(number of processor may needed; int), sorted_chr_values(sizes of peak files, sorted based on chr number (decending); list), sorted_chr_keys (chromosome number, decending sorted; list), chr_count (number of chromosomes need to run; int), cur_wig_subfolder (sample type; list), path_to_wig_files(/path/to/wig/file/; str), output_dir (/path/to/out/sample/; str), num_windows (number of windows want to run; int), dist_each_dir (distance to each directions; int)
#output: 
#e.g.: ProcessorSchedule(num_processors, sorted_chr_values, sorted_chr_keys,chr_count,cur_wig_subfolder, path_to_wig_files, output_dir, num_windows, dist_each_dir)
def ProcessorSchedule(num_processors,sorted_chr_values, sorted_chr_keys, chr_count, cur_wig_subfolder, path_to_wig_files, output_dir, num_windows, dist_each_dir):
    sizes = {}
    processors = {}
    children = []
    for i in range(1, num_processors+1):
        sizes[i] = sorted_chr_values[0]
        processors.setdefault(i,[]).append(sorted_chr_keys[0])
        sorted_chr_keys.pop(0)
        sorted_chr_values.pop(0)

    max_size = chr_count - num_processors #get the number of chromomes over the number of processors available 
    counter = 1
    while counter <= max_size:
        #put sorted keys (chr number) in a list (ks)
        ks = list(dict(sorted(sizes.items(), key=operator.itemgetter(0))).keys())
        processors.setdefault(ks[counter-1],[]).append(sorted_chr_keys[0])
        if sizes.get(ks[counter-1]): #in case the is key not in sizes before
            sizes[ks[counter-1]] = sorted_chr_values[0]
        else:
            sizes[ks[counter-1]] = sizes[ks[counter-1]] + sorted_chr_values[0] #if exist key, add tow values
        sorted_chr_keys.pop(0)
        sorted_chr_values.pop(0)
        counter+=1

    chr_sets = [] #after for loop finished, it should be a list of list
    for chr_i in range(max_size+1):
        temp_l = []
        for count in range(1, num_processors+1):
            #if chr_i greater than number of values of count  
            try:
                processors[count][chr_i]
            except NameError:
                next
            else:
                temp_l.append(processors[count][chr_i])
        chr_sets.append(temp_l)

    #start program
    for i in range(len(chr_sets)):
        chrs = chr_sets[i]
        for count in range(1, num_processors+1):
            try:
                pid = os.fork()
            except OSError:
                sys.stderr.write("Couldn't fork, no child process create\n")
                continue

            if pid == 0:
                #child 
                print("In the child process {} that has the PID {}".format(count, os.getpid()))
            
                if chrs.get(count-1):
                    full_path_to_sample_wigs = path_to_wig_files+cur_wig_subfolder+"/"
                    args = [count, full_path_to_sample_wigs, output_dir, num_windows, dist_each_dir, chrs[count-1]] #array of chroms assigned to pid, pid, cur_wig_subfolder, path_to_wig_files, output directory
                
                    runWindowMedians(args)
                    exit()
            else:
                #parent 
                children.append(pid)
        for c in children:
            tmp = os.waitpid(c, 0)
    return args





