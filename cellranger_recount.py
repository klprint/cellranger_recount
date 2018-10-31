###
# Created: 26.06.2018
# By: Kevin Leiss
###

import pandas as pd
import numpy as np
import sys
import re
from subprocess import check_output
from multiprocessing.dummy import Pool as ThreadPool
import itertools
import HTSeq

# create a GTF dict
def read_gtf(filepath, subset_kind = "transcript"):
    '''Generates a dictionary which contains the gene_ids as keys.
    Stored informations are: Chromosome, type, start-position, stop-position and strandness'''
    gtf_file = HTSeq.GFF_Reader(filepath)

    subset = HTSeq.GenomicArrayOfSets("auto", stranded = True)

    for feature in gtf_file:
        if feature.type == subset_kind:
            subset[feature.iv] = feature.name

    return(subset)


def process_alignment(alnmt, featureset):
    '''Processes HTSeq alignments and uses a parsed featureset (generate by read_gtf)
    to return the Barcode, UMI and the aligning gene.
    Returns None if there was no gene annotated'''
    if alnmt.aligned:
            if alnmt.aQual == 255:
                optional_fields = [keys[0] for keys in alnmt.optional_fields]
                if "MM" not in optional_fields and "CB" in optional_fields and "UB" in optional_fields:
                    features = list(featureset[alnmt.iv].steps())
                    features = [feature for feature in features if feature[1] != set()]

                    if len(features) > 0:

                        return(alnmt.optional_field("CB"), alnmt.optional_field("UB"), tuple(feature[1] for feature in features))



def get_aligned_molecules(bamfilepath, featureset, umi = True):
    '''Returns a array of: 
      1. list of tuples of Barcode, UMI, aligned gene
      2. int of alignments which did not pass filters'''
    bamfile = HTSeq.BAM_Reader(bamfilepath)

    molecules_out = []
    umi_out = set()
    not_passed = 0

    for alnmt in bamfile:
        if alnmt.aligned:
            mol_out = process_alignment(alnmt, featureset)
            if mol_out == None:
                not_passed += 1
                continue
            if umi:
                umi_out.add(mol_out)
            else:
                molecules_out.append(mol_out)

    if umi:
        return [u for u in umi_out], not_passed
    else:
        return molecules_out, not_passed


def parse_molecules(aligned_molecules):
    outdict = dict()
    for aln in aligned_molecules:
        genes = aln[2]

        if aln[0] not in outdict.keys():
            outdict[aln[0]] = dict()
            outdict[aln[0]] = {key: 1 for key in genes}
        else:
            for gene in genes:
                if gene not in outdict[aln[0]].keys():
                    outdict[aln[0]][gene] = 1
                else:
                    outdict[aln[0]][gene] += 1
    return outdict


def make_sparse_mtx(parsed_molecules, outfolder):
    barcodes = [key for key in parsed_molecules.keys()]
    genes = set()
    countsum = 0

    # for key, entry in parsed_molecules.items():
    #     reg_genes = [gene for gene in entry.keys()]
    #     for gene in reg_genes:
    #         genes.add(gene)

    print("Parsing coodrinates")
    for key, entry in parsed_molecules.items():
        for gene, count in entry.items():
            genes.add(gene)
            countsum = countsum + count

    genes = list(genes)

    print("Writing matrix to " + outfolder)
    with open(outfolder + "/matrix.mtx", "a") as out_mtx:
        out_mtx.write("%%MatrixMarket matrix coordinate integer general\n")
        out_mtx.write("%\n")
        out_mtx.write(str(len(genes)) + " " + str(len(barcodes)) + " " + str(countsum) + "\n")


        for bc, entry in parsed_molecules.items():
            bc_coord = barcodes.index(bc) +1  # Make it 1-based counting

            for gene, count in parsed_molecules[bc].items():
                gene_coord = genes.index(gene) +1  # Make it 1-based
                out_mtx.write(str(gene_coord) + " " + str(bc_coord) + " " + str(count) + "\n")

    print("Writing genes.tsv")
    with open(outfolder + "/genes.tsv", "a") as out_genes:
        for gene in genes:
            out_genes.write(gene + "\t" + gene + "\n")

    print("Writing barcodes.tsv")
    with open(outfolder + "/barcodes.tsv", "a") as out_barcodes:
        for bc in barcodes:
            out_barcodes.write(bc + "\n")








if __name__ == "__main__":
    pass
