import matplotlib as mpl
import argparse
mpl.use("Agg")
import matplotlib.pyplot as plt
import multiprocessing as mp
from multiprocessing import pool
import sys, re, os
import numpy as np 
import bpRNA_CosMoS_module as cosmos

def multi_process(arg_list, run_function, n_cores):
    #print(arg_list)
    pool = mp.Pool(n_cores)  # create a pool of workers w/ n cores
    results = pool.map(run_function, arg_list)  # use pool w/ map to run affine samplers in parallel
    pool.close()
    return results

def get_struct_info(st_file):
    # Get structure array and db information
    with open(st_file, 'r') as st_info:
        count = 0
        for i, line in enumerate(st_info):
            if line[0] != "#":
                count += 1
            if count == 2:
                db = line.strip()
            if count == 3:
                ss = line.strip()
    return ss, db

#Generate updated ss array
def edit_ss_array(ss, dotbracket):
    ss_edit = ''
    right_symbols = [")", "]", ">", "}"]
    left_symbols = ["(", "[", "<", "{"]
    for i, s in enumerate(ss):
        if (dotbracket[i] in right_symbols or dotbracket[i].islower()) and s == "S":
            ss_edit += "R"
        elif (dotbracket[i] in left_symbols or dotbracket[i].isupper()) and s == "S":
            ss_edit += "L"
        else:
            ss_edit += s
    return ss_edit

def run_get_similarity_score(arg_list):
    score_results = []
    for arg in arg_list:
        results = get_similarity_score(arg)
        score_results.append(results)
    return score_results 

def get_similarity_score(arg_list):
    name_1, name_2, ss_1, ss_2 = arg_list
    name_combination = name_1 + "\t" + name_2
    if kmer_distances:
        sim_score = cosmos.get_similarity_score(ss_1, ss_2, k, kmer_distances) 
    else:
        sim_score = cosmos.get_similarity_score(ss_1, ss_2, k) 
    return (sim_score, name_combination)

def get_paths(file_paths):
    paths = []
    with open(file_paths, "r") as paths_data:
        for line in paths_data:
            paths.append(line.strip())
    return paths

def get_structure_dict(paths):
    strct_dict = {}
    for RNA_file in paths:
        name = RNA_file.split("/")[-1].split(".")[0]
        file_type = RNA_file.split("/")[-1].split(".")[-1]
        if file_type != "st":
            if file_type == "dbn" or file_type == "DBN":
                os.system("perl bpRNA.pl " + RNA_file)
                RNA_file = name + ".st"
            else:
                print("Incorrect file type used, try dbn or st instead")
                sys.exit()
        ss, dbn = get_struct_info(RNA_file)
        ss = edit_ss_array(ss, dbn)
        strct_dict[name] = ss
    return strct_dict

def get_arg_list(strct_dict):
    arg_list = []
    name_combinations = []
    for name_1 in strct_dict:
        ss_1 = strct_dict[name_1]
        for name_2 in strct_dict:
            if name_1 + "\t" + name_2 not in name_combinations and name_2 + "\t" + name_1 not in name_combinations and name_1 != name_2:
                ss_2 = strct_dict[name_2]
                name_combinations = np.append(name_combinations, name_1 + "\t" + name_2)
                arg_list.append((name_1, name_2, ss_1, ss_2)) 
    return arg_list

def write_files(output_file_name, score_results):
    with open(output_file_name, 'w') as similarity_file:
        similarity_file.write("name_1\tname_2\tsimilarity_score")
        for data in score_results:
            sim_score, name_combination = data
            similarity_file.write("\n" + name_combination + "\t" + str(sim_score))

########
##MAIN##
########

# Add argument flags and required inputs
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--path_list_file", required=True, help="A text file containing a list of .st files or .dbn file paths", type=str)
parser.add_argument("-a", "--fuzzy", required=False, help="Use fuzzy option (True, False)", choices=[True, False], default=False, type=bool)
parser.add_argument("-o", "--output_base_name", required = False, default = "test", help="Output file base name containing alignment score results", type=str)
parser.add_argument("-c", "--cores", required = False, default = 10, help="Cores inputs the number of cores to use when using multiprocessing", type=int)
parser.add_argument("-p", "--mp", required = False, help="option to run multiprocessing", choices=[True, False], default=False, type=bool)
args = parser.parse_args()

# Asign variables from input arguments
file_paths = args.path_list_file
use_fuzzy = args.fuzzy
name = args.output_base_name
n_cores = args.cores
multiprocessing = args.mp
# Function calls
if use_fuzzy:
    kmer_edit_dist_file = "kmer_pc_dict_max_dist_2.npy"
    kmer_distances = np.load(kmer_edit_dist_file, allow_pickle='TRUE').item()
    k = 9
else:
    kmer_distances = False
    k = 10

#name = "bpRNA_align_benchmark_data"
paths = get_paths(file_paths)
output_file_name = name + "_bpRNA_CosMoS_scores.txt"
strct_dict = get_structure_dict(paths)
arg_list = get_arg_list(strct_dict)
if multiprocessing:
    score_results = multi_process(arg_list, get_similarity_score, n_cores)  
else:
    score_results = run_get_similarity_score(arg_list)
write_files(output_file_name, score_results)
