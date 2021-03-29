#!/usr/bin/env python

import sys
import subprocess
from collections import defaultdict

#####
# This script calculates coverage in sliding windows from BAM(s) and outputs compiled information from all BAMs input.
# Note: "position" in output refers to the beginning of the window.
#####

if len(sys.argv) < 8:
    print("Usage: slidingWindowCoverage.py path_to_samtools reference_header window_size step_size start end <bam1> <bam2> ... <bamN>")
    sys.exit(0)

sam_path = sys.argv[1]
ref = sys.argv[2]
w = int(sys.argv[3])
s = int(sys.argv[4])
start = int(sys.argv[5])
end = int(sys.argv[6])
bams = sys.argv[7:]


# make bed file of regions to calculate coverage
print("Making bed file of windows to calculate coverage ...")
bed = "slidingWindowCoverage.bed"
with open(bed,"w") as out:
    for x in range(start,end,s):
        out.write("\t".join([ref,str(x),str(x+w),"cov"])+"\n")

# run samtools bedcov on each bam file
print("Running samtools bedcov ...")
outs = []
for bam in bams:
    out = bam.split(".")[0]+"_cov.txt"
    outs.append(out)
    with open(out,"w") as f:
        subprocess.call([sam_path,"bedcov",bed,bam], stdout=f)

# combine data from multiple BAMs
print("Combining coverage data from bams ...")
covDict = defaultdict(list)
for data in outs:
    with open(data,"r") as f:
        for line in f:
            info = line.strip().split("\t")
            pos = info[1]
            try:
                cov = int(info[4])
            except ValueError:
                print(line)
            cov = round(cov/w,0)
            covDict[pos].append(cov)

# write output
print("Writing to output ...")
with open("slidingWindowCoverage_"+str(w)+"w"+str(s)+"s"+".txt","w") as out:
    header = ["position"]+[x.split("_")[0] for x in outs]
    out.write("\t".join(header)+"\n")
    for key,value in covDict.items():
        out.write(key+"\t"+"\t".join([str(x) for x in value])+"\n")
