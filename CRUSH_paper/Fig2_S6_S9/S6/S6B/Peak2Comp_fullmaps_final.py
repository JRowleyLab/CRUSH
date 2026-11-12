#!/usr/bin/env python
from scipy import stats as scstat
import numpy as np
import os
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
import argparse, argcomplete
import subprocess
import pyBigWig
from scipy.stats import mstats
from scipy.stats.mstats import gmean as gmean
import warnings
import pandas as pd

parser = argparse.ArgumentParser(description='Tool to correlate compartments with bigwigs\n\n')

parser.add_argument("-e", dest="efile", help="Input file of eigenvector in 4th column\n\n", required=True)
parser.add_argument('-b', dest='bfile', help="Input bigwig file\n\n", required=True)
parser.add_argument("-o", dest="ofile", help="path to output files\n\n", metavar="ofile", required=True)


argcomplete.autocomplete(parser)
args = parser.parse_args()


bfile=str(args.bfile)
efile=str(args.efile)
ofile=str(args.ofile)
chrlist=[]
startlist=[]
endlist=[]
compvaluelist=[]
with open(efile, "r") as e:
    for lines in e:
        lines = lines.strip('\n')
        line=lines.split("\t")
        chrlist.append(str(line[0]))
        startlist.append(int(line[1]))
        endlist.append(int(line[2]))
        compvaluelist.append(float(line[3]))

bw = pyBigWig.open(bfile)
geoms = []
means = []
medians = []
for comp in tqdm(range(0,len(chrlist))):
    mychr = "chr" + str(chrlist[comp])

    rlist=bw.values(mychr,startlist[comp],endlist[comp])
    rlist2=[]
    for item in rlist:
        if item == 0:
            continue
        else:
            rlist2.append(item)
    if len(rlist2) > 0:
        warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')
        warnings.filterwarnings(action='ignore', message='Mean of empty slice')
        geoms.append(gmean(rlist2))
        means.append(np.nanmean(rlist2))
        medians.append(np.nanmedian(rlist2))
    else:
        geoms.append(0)
        means.append(0)
        medians.append(0)
geoms2 = np.nan_to_num(geoms)
means2 = np.nan_to_num(means)
medians2 = np.nan_to_num(medians)

data = {"chr":chrlist, "start":startlist, "end":endlist, "values":compvaluelist, "means1":means, "medians1":medians, "geoms1":geoms, "means2":means2, "medians2":medians2, "geoms2":geoms2}

df = pd.DataFrame(data)

df2= df[df['means1'].notna()]


df2.to_csv(ofile, index=False, header=True, sep="\t")
