#!/usr/bin/env python

import argparse
import os.path
from collections import defaultdict
from datetime import datetime

def get_args():
    """
    Handle the command line arguments.
    """
    parser = argparse.ArgumentParser(description="Converts Popoolation2 sync file to TSV with allele frequencies. \
        Filters mutations based on bed file,\
        allele count, coverage and frequency across samples. NOTE: only considers ACTG SNPs (no N's, no deletions).")
    parser.add_argument("syncfile",help="sync file from Popoolation2")
    parser.add_argument("--bed",help="bed file specifying regions to remove from analysis")
    parser.add_argument("--min-count",default=5,type=int,help="minimum count of alternate allele. default= 5")
    parser.add_argument("--min-coverage",default=30,type=int,help="minimum coverage (across all samples) for SNP to be \
        called. default= 30")
    parser.add_argument("--min-freq",default=5,type=int,help="minimum variant frequency to be considered - any variant \
        below this frequency will be changed to zero. default= 5%%")
    parser.add_argument("--output",default="popoolation_sync.tsv",help="Filename to write output to. default= 'popoolation_sync.tsv'")
    return parser.parse_args()

def getBedCoor(bedfile):
    """
    Reads bed file into list of coordinates.
    """
    bedCoor = []
    with open(bedfile, "r") as f:
        for line in f:
            info = line.strip().split("\t")
            for i in range(int(info[1]), int(info[2])+1):
                bedCoor.append(i)
    return(bedCoor)

def alleleFreq(syncfile, filter_coord, outfile, mincount, mincov, minfreq):
    """
    Reads sync file, filters mutations based on position, allele count, coverage and variant frequency.
    """
    print("Start: "+str(datetime.now()))
    freq_dict = defaultdict(lambda: defaultdict(list))
    base_dict = {0:"A",1:"T",2:"C",3:"G",4:"DEL"}
    rev_base_dict = {"A":0,"T":1,"C":2,"G":3,"DEL":4}
    with open(syncfile,"r") as f:
        print("Parsing sync file ...")
        for line in f:
            info = line.strip().split("\t")
            pos = info[1]
            if int(pos) % 100000 == 0:
                print("Processed "+pos+" positions")
            ref = info[2]
            ref_index = rev_base_dict[ref]
            time_points = info[3:]
            num_time_points = len(time_points)
            # check if all time points meet minimum coverage requirements & position not in coordinates to filter out
                # checking that sum of allele counts is >= min-coverage for all sequenced time points
            if (all(sum([int(x) for x in y.split(":")])>= mincov for y in time_points)) and (pos not in filter_coord):
                pos_key = pos+"_"+ref
                for point in time_points:
                    total = sum([int(x) for x in point.split(":")])
                    base_counts = point.split(":")
                    # only looking at counts of ATCG - no N's, no deletions
                    base_index = -1
                    positions = [0,1,2,3,5]
                    for x in positions:
                        base = base_counts[x]
                        base_index += 1
                        # check if alternate allele counts meet minimum count requirements
                        if int(base) >= mincount:
                            freq = round((int(base)/total*100),0)
                            # check if variant frequency is above minimum threshold
                            if freq < minfreq:
                                freq = 0      
                        else:
                            freq = 0
                        alt = base_dict[base_index]
                        freq_dict[pos_key][alt].append(freq)
    # header for output file
    header = "\t".join(["pos","ref","alt"]+["t"+str(n) for n in range(1,num_time_points+1)])
    with open(outfile, "w") as out:
        print("Filtering and writing to output ...")
        out.write(header+"\n")
        for key,value in freq_dict.items():
            pos = key.split("_")[0]
            ref = key.split("_")[1]
            alt = []
            time_dict = defaultdict(list)
            for x,y in value.items():
                y = [int(i) for i in y]
                # exclude frequency of reference allele & exclude alleles with 0% freq
                if (not x == ref) and (not all(i == 0 for i in y)):
                    alt.append(x)
                    pos_count = 0
                    for i in y:
                        pos_count += 1
                        time_dict[pos_count].append(i)
            if len(alt) != 0:
                alt = ",".join(alt)
                out.write("\t".join([pos,ref,alt]))
                for time,freqs in time_dict.items():
                    out.write("\t"+",".join([str(x) for x in freqs]))
                out.write("\n")
    print("End: "+str(datetime.now()))



args = get_args()

if args.bed is not None:
    bad_coor = getBedCoor(args.bed)
else:
    bad_coor = []

alleleFreq(args.syncfile, bad_coor, args.output, args.min_count ,args.min_coverage ,args.min_freq)

