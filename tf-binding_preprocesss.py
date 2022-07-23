import sys
import pandas as pd
from collections import defaultdict
from pandas_plink import read_plink1_bin
from argparse import ArgumentParser
import timeit
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning
import json

simplefilter("ignore", category=ConvergenceWarning)
# general params
CUTOFF_PVALUE = 0.05
NUM_SHUFFLES = 100
NUM_ROBUST_ITER = 5
FRACTION_ROBUST_SAMPLES = 0.9

# read the names of input and output files
parser = ArgumentParser(description='Get parameters')
parser.add_argument('-s', dest='snpSiftFile', type=str,
                   help='SNPs after filtered by SnpSift')
parser.add_argument('-g', dest='genotypeFile', type=str,
                   help='Genotype file')
parser.add_argument('-o', dest='oVCF', type=str,
                   help='Output VCF file')
parser.add_argument('-otf', dest='TF2pos', type=str,
                   help='Output TF2pos file')

if len(sys.argv) <9:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

start = timeit.default_timer()

# Mapping from gene to its TFs
tfSet = set()
tf2Genes = defaultdict(set)
with open('TF2Gene.txt', 'r') as f:
    line = f.readline()  # header
    for line in f:
        cols = line.strip().split('\t')
        TF = cols[0]
        gene = cols[1]
        tfSet.add(TF)
        tf2Genes[TF].add(gene)


# Mapping between gene codes and gene name
pos2TF = {}
pos2rsid = {}
with open(args.snpSiftFile, 'r') as f:
    line = f.readline()
    while line.startswith('#'):
        line = f.readline()
    for i, line in enumerate(f):
        cols = line.strip().split('\t')
        chrom = "chr" + str(cols[0])
        pos = str(cols[1])
        rsid = str(cols[2])
        info = cols[-1].split(';')
        for field in info:
            if field.startswith('GENEINFO'):
                gene = field.split('=')[1].split(':')[0]
                if gene in tfSet:
                    pos2TF[chrom + ":" + pos] = gene
                    pos2rsid[chrom + ":" + pos] = rsid

outVCF = open(args.oVCF, 'w')
outTf2Pos = open(args.TF2pos, 'w')
outTf2Pos.write('TF\tposition\n')
with open(args.genotypeFile, 'r') as f:
    line = f.readline()
    index = 0
    while line.startswith('#'):
        outVCF.write(line)
        line = f.readline()
    for i, line in enumerate(f):
        cols = line.strip().split('\t')
        chrom = str(cols[0])
        pos = str(cols[1])
        full_pos = chrom + ":" + pos
        if full_pos in pos2TF:
            outVCF.write(line)
            outTf2Pos.write(f'{pos2TF[full_pos]}\t{pos2rsid[full_pos]}\t{str(index)}\n')
            index = index+1

end = timeit.default_timer()
print('Finished mapping genes to TFs. Took ' + str(int(end-start)) + ' seconds')
