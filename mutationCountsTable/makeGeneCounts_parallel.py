#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
<<<<<<< HEAD
# script: makeGeneCounts_parallel.py
# author: Lincoln 
=======
# script: makeGeneCounts.py
# author: Lincoln
>>>>>>> b6191860c404154094b8aa86beae91fa7578464d
# date: 12.18.18
#
# Creates a gene/cell table for the mutations found in a given
# population of cells. Trying to implement parallelization here
# Run this on a big-ass machine
#
# usage:
#			python3 makeGeneCounts_parallel.py
#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
import argparse
import numpy as np
import VCF # comes from Kamil Slowikowski
import os
import os.path
import pickle
import csv
import pandas as pd
import re
import sys
import multiprocessing as mp
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def getFileNames(input_dir):
        # Get file names based on the specified path
        files = []
        for file in os.listdir(input_dir):
                if file.endswith(".vcf"):
                        fullPath = (os.path.join(input_dir, file))
                        files.append(fullPath)

        return files


def getGenomePos(sample):
        # Returns a genome position sting that will match against the ones w/in COSMIC db
        chr = str(sample[0])
        chr = chr.replace("chr", "")
        pos = int(sample[1])
        ref = str(sample[3])
        alt = str(sample[4])

        if (len(ref) == 1) and (len(alt) == 1): # most basic case
                secondPos = pos
        elif (len(ref) > 1) and (len(alt) == 1):
                secondPos = pos + len(ref)
        elif (len(alt) > 1) and (len(ref) == 1):
                secondPos = pos + len(alt)
        else: # multibase-for-multibase substitution
                alts = alt.split(",")
                assert len(alts) > 1
                len_alts = [len(alt) for alt in alts]
                secondPos = pos + max(len_alts)

        genomePos = chr + ':' + str(pos) + '-' + str(secondPos)
        return genomePos


def getGeneName(posString):
        # want to return the gene name from a given genome position string
        # (ie. '1:21890111-21890111'), by querying the hg38-plus.gtf
        global m19_lookup

        # work on posString
        chrom = posString.split(':')[0]
        posString_remove = posString.split(':')[1]
        lPosition = int(posString_remove.split('-')[0])
        rPosition = int(posString_remove.split('-')[1])

        # work on m19_gtf
        chromStr = 'chr' + str(chrom)
        gene_names = m19_lookup[chromStr].overlap(rPosition, lPosition+1)
        gene_names = list({gname.data for gname in gene_names})
        if len(gene_names) > 0:
                return gene_names
        else:
                return ""


def getGeneCellMutCounts(f):
        # Creates dictionry obj where every key is a cell and every value is
        # a list of the genes we found mutations in for that cell.
        tup = [] # not really a tuple, just a list, i guess

        cell = os.path.basename(f)
        cell = cell.replace(".vcf", "")
        print(cell) # to see where we are

        df = VCF.dataframe(f)
        genomePos_query = df.apply(getGenomePos, axis=1) # apply function for every row in df

        shared = list(set(genomePos_query)) # genomePos_query (potentially) has dups
        sharedGeneNames = [f for e in shared for f in getGeneName(e)]
        tup = [cell, sharedGeneNames]

        return(tup)


def formatDataFrame(raw_df):
        # logic for creating the cell/mutation counts table from the raw
        # output that getGeneCellMutCounts provides
        cellNames = list(raw_df.index)

        genesList = []
        for i in range(0, raw_df.shape[0]):
                currList = list(raw_df.iloc[i].unique()) # unique genes for curr_cell

                for elm in currList:
                        if elm not in genesList:
                                genesList.append(elm)

        genesList1 = pd.Series(genesList)

        df = pd.DataFrame(columns=genesList1, index=cellNames) # initialize blank dataframe
        for col in df.columns: # set everybody to zero
                df[col] = 0

        for i in range(0,raw_df.shape[0]):
                currCell = raw_df.index[i]
                currRow = raw_df.iloc[i]

                for currGene in currRow:
                        df[currGene][currCell] += 1

        return(df)


def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument("--ref-gtf",
                            type=argparse.FileType('r'),
                            required=True)
        parser.add_argument("--input-dir",
                            required=True)
        parser.add_argument("--output",
                            type=argparse.FileType('w'),
                            default="geneCellMutationCounts_test_toggle.csv")
        parser.add_argument("--nprocs", type=int, default=16)
        parser.add_argument("--interval-tree",
                            type=argparse.FileType("r"),
                            required=True)

        return parser.parse_args()


def main():
        global m19_lookup

        args = parse_args()

        m19_lookup = pickle.load(open(args.interval_tree.name, 'rb'))
        fNames = getFileNames(args.input_dir)

        print('creating pool')

        p = mp.Pool(processes=args.nprocs)

        try:
                cells_list = p.map(getGeneCellMutCounts, fNames, chunksize=1) # default chunksize=1
        finally:
                p.close()
                p.join()

        cells_dict = {}

        for item in cells_list:
            cells_dict.update({item[0]:item[1]})

        print('writing file')

        filterDict_pd = pd.DataFrame.from_dict(cells_dict, orient="index") # orient refers to row/col orientation
        filterDict_format = formatDataFrame(filterDict_pd)
        filterDict_format.to_csv(args.output.name)


if __name__ == "__main__":
        main()
