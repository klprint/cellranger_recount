#!/usr/bin/env python
###
# Created: 26.06.2018
# By: Kevin Leiss
###

import pandas as pd
import numpy as np
import sys
import time
import re
import itertools
from subprocess import check_output
from multiprocessing import Pool as ThreadPool
import gzip

# create a GTF dict
def read_gtf(filepath, subset_kind = "transcript"):
    '''Generates a dictionary which contains the gene_ids as keys.
    Stored informations are: Chromosome, type, start-position, stop-position and strandness'''

    if filepath.endswith(".gz"):
        with gzip.open(filepath, "r") as gtfile:
            gtfdict = dict()
            for line in gtfile:
                line = line.decode("utf-8")
                if(line.startswith("#")):
                    continue
                else:
                    entries = line.split("\t")
                    
                    chrom = entries[0]
                    kind = entries[2]

                    if(kind != subset_kind):
                        continue

                    fstart = entries[3]
                    fstop = entries[4]
                    strand = entries[6]
                    meta = entries[8].split("; ")

                    metadict = {m.split(" ")[0]: m.split(" ")[1].replace("\"", "") for m in meta}

                    if(metadict["gene_id"] in gtfdict.keys()):
                        gtfdict[metadict["gene_id"]].append((chrom, fstart, fstop, strand))
                    else:
                        gtfdict[metadict["gene_id"]] = [(chrom, fstart, fstop, strand)]
    else:
        with open(filepath, "r") as gtfile:
            gtfdict = dict()
            for line in gtfile:
                if(line.startswith("#")):
                    continue
                else:
                    entries = line.split("\t")
                    
                    chrom = entries[0]
                    kind = entries[2]

                    if(kind != subset_kind):
                        continue

                    fstart = entries[3]
                    fstop = entries[4]
                    strand = entries[6]
                    meta = entries[8].split("; ")

                    metadict = {m.split(" ")[0]: m.split(" ")[1].replace("\"", "") for m in meta}

                    if(metadict["gene_id"] in gtfdict.keys()):
                        gtfdict[metadict["gene_id"]].append((chrom, fstart, fstop, strand))
                    else:
                        gtfdict[metadict["gene_id"]] = [(chrom, fstart, fstop, strand)]

    return gtfdict

                

def call_samtools(bamfilepath, chrom, start, stop, strand, alignment_qual="255", strandness="sense",):
    # if strandness != "sense" or strandness != "antisense":
    #     raise ValueError("strandness needs to be either sense or antisense.")

    if strandness == "antisense":
        if strand == "+":
            strand = "-"
        else:
            strand = "+"

    if strand == "+":
        aligning_reads = check_output(["samtools", "view", "-q", alignment_qual, "-F 16", bamfilepath, chrom + ":" + start + "-" + stop]).decode("utf-8").split("\n")
    else:
        aligning_reads = check_output(["samtools", "view", "-q", alignment_qual, "-f 16", bamfilepath, chrom + ":" + start + "-" + stop]).decode("utf-8").split("\n")

    del aligning_reads[-1]
    return aligning_reads



def get_read_metadata(readstring):
    read_info = readstring.split("\t")
    meta_block = read_info[11:len(read_info)]
    meta_block = {re.split(":.:",m)[0]: re.split(":.:",m)[1] for m in meta_block}

    return meta_block

def filter_reads(readlist, cr_altered = False, unique = True):
    out_readslist = list()

    for read in readlist:
        meta_block = get_read_metadata(read)

        if(meta_block["NH"] == "1" and "MM" not in meta_block.keys() and "UB" in meta_block.keys() and "CB" in meta_block.keys()):
            out_readslist.append(read)

    return(out_readslist)

def assign_read_identitiy(read):
    meta_block = get_read_metadata(read)
    return (meta_block["CB"], meta_block["UB"])


def call_samtools_by_gtfdict(bamfilepath, gene, gtfdict, alignment_qual="255", strandness="sense", cr_altered = False, unique = True):
    entries = gtfdict[gene]
    all_reads = set()

    for entry in entries:
        chrom, start, stop, strand = entry
        found_reads = call_samtools(bamfilepath, chrom, start, stop, strand, alignment_qual=alignment_qual, strandness=strandness)
        if(len(found_reads) > 0):
            for read in found_reads:
                all_reads.add(read)


    filtered_reads = filter_reads(list(all_reads), cr_altered = cr_altered, unique = unique)
    id_reads = [assign_read_identitiy(read) for read in filtered_reads]

    id_reads = [(gene, cell, umi) for cell, umi in id_reads]
    
    return id_reads


def updt(total, progress):
    """
    Displays or updates a console progress bar.

    Original source: https://stackoverflow.com/a/15860757/1391441
    """
    barLength, status = 20, ""
    progress = float(progress) / float(total)
    if progress >= 1.:
        progress, status = 1, "\r\n"
    block = int(round(barLength * progress))
    text = "\r[{}] {:.0f}% {}".format(
        "#" * block + "-" * (barLength - block), round(progress * 100, 0),
        status)
    sys.stdout.write(text)
    sys.stdout.flush()


def check_reads_per_feature(bamfilepath, gtf_dict, alignment_qual="255", strandness="sense", ncores = 4, umi = True):
    genes = list(gtf_dict.keys())

    out = list()

    with ThreadPool(ncores) as pool:
        for o in pool.starmap(call_samtools_by_gtfdict, zip(itertools.repeat(bamfilepath), genes, itertools.repeat(gtf_dict), itertools.repeat(alignment_qual), itertools.repeat(strandness))):
            out = out + o

    if(umi):
        out = list(set(out))
    return out
    

def generate_sparse_matrix(gene_bc_umi_tuplelist, outdir):
    genes = set()
    barcodes = set()

    n_entries = len(gene_bc_umi_tuplelist)

    i=1
    for entry in gene_bc_umi_tuplelist:
        genes.add(entry[0])
        barcodes.add(entry[1])
        updt(n_entries, i)
        i+=1


# def process_alignment(alnmt, featureset):
#     '''Processes HTSeq alignments and uses a parsed featureset (generate by read_gtf)
#     to return the Barcode, UMI and the aligning gene.
#     Returns None if there was no gene annotated'''
#     if alnmt.aligned:
#             if alnmt.aQual == 255:
#                 optional_fields = [keys[0] for keys in alnmt.optional_fields]
#                 if "MM" not in optional_fields and "CB" in optional_fields and "UB" in optional_fields:
#                     features = list(featureset[alnmt.iv].steps())
#                     features = [feature for feature in features if feature[1] != set()]

#                     if len(features) > 0:

#                         return(alnmt.optional_field("CB"), alnmt.optional_field("UB"), tuple(feature[1] for feature in features))



# def get_aligned_molecules(bamfilepath, featureset, reads = False):
#     '''Returns a array of: 
#       1. list of tuples of Barcode, UMI, aligned gene
#       2. int of alignments which did not pass filters'''
#     bamfile = HTSeq.BAM_Reader(bamfilepath)

#     molecules_out = []
#     umi_out = set()
#     not_passed = 0

#     n_processed = 1
#     for alnmt in bamfile:
#         if alnmt.aligned:
#             mol_out = process_alignment(alnmt, featureset)
#             if mol_out == None:
#                 not_passed += 1
#                 continue
#             if reads:
#                 molecules_out.append(mol_out)
#             else:
#                 umi_out.add(mol_out)
                
#         if n_processed % 500000 == 0:
#             print(str(n_processed) + " alignments processed")
#         n_processed += 1

#     if reads:
#         return molecules_out, not_passed
#     else:
#         return [u for u in umi_out], not_passed


# def parse_molecules(aligned_molecules):
#     outdict = dict()
#     for aln in aligned_molecules:
#         genes = aln[2]

#         if aln[0] not in outdict.keys():
#             outdict[aln[0]] = dict()
#             outdict[aln[0]] = {key: 1 for key in genes}
#         else:
#             for gene in genes:
#                 if gene not in outdict[aln[0]].keys():
#                     outdict[aln[0]][gene] = 1
#                 else:
#                     outdict[aln[0]][gene] += 1
#     return outdict


# def make_sparse_mtx(parsed_molecules, outfolder):
#     barcodes = [key for key in parsed_molecules.keys()]
#     genes = set()
#     countsum = 0

#     # for key, entry in parsed_molecules.items():
#     #     reg_genes = [gene for gene in entry.keys()]
#     #     for gene in reg_genes:
#     #         genes.add(gene)

#     print("Parsing coodrinates")
#     for key, entry in parsed_molecules.items():
#         for gene, count in entry.items():
#             genes.add(gene)
#             countsum = countsum + count

#     genes = list(genes)

#     print("Writing matrix to " + outfolder)
#     with open(outfolder + "/matrix.mtx", "a") as out_mtx:
#         out_mtx.write("%%MatrixMarket matrix coordinate integer general\n")
#         out_mtx.write("%\n")
#         out_mtx.write(str(len(genes)) + " " + str(len(barcodes)) + " " + str(countsum) + "\n")


#         for bc, entry in parsed_molecules.items():
#             bc_coord = barcodes.index(bc) +1  # Make it 1-based counting

#             for gene, count in parsed_molecules[bc].items():
#                 gene_coord = genes.index(gene) +1  # Make it 1-based
#                 out_mtx.write(str(gene_coord) + " " + str(bc_coord) + " " + str(count) + "\n")

#     print("Writing genes.tsv")
#     with open(outfolder + "/genes.tsv", "a") as out_genes:
#         for gene in genes:
#             out_genes.write(gene + "\t" + gene + "\n")

#     print("Writing barcodes.tsv")
#     with open(outfolder + "/barcodes.tsv", "a") as out_barcodes:
#         for bc in barcodes:
#             out_barcodes.write(bc + "\n")








if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Recounting cellranger BAM files')
    parser.add_argument("bam", 
        metavar="BAM", 
        type=str, 
        default = None,
        help='Path to the to be analyzed BAM file')

    parser.add_argument("gtf",
        metavar="GTF",
        type=str,
        default=None,
        help="Path to the genome annotation (GTF/GFF-format)")

    parser.add_argument("outdir", 
        metavar="OUTDIR", 
        type = str, 
        default = None,
        help = "Directory where the output matrices should be stored")


    parser.add_argument("-f", "--feature",
        metavar="FEATURE",
        type=str,
        default="transcript",
        help="What kind of feature should be counted? default: transcript")


    parser.add_argument("--reads", 
        action = "store_true",
        help = "If set, the returned matrix will be readcounts instead of UMI counts")


    args = parser.parse_args()

    #print(args)


    import datetime as dt
    import random
    import os

    logfile = "log_" + str(random.randint(0,1000)) + ".txt"
    sys.stdout = open(logfile, 'w',  buffering = 1)

    print("StdOut will be redirected into file: " + logfile)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    else:
        if os.path.exists(args.outdir + "/matrix.mtx") or os.path.exists(args.outdir + "/genes.tsv") or os.path.exists(args.outdir + "/barcodes.tsv"):
            raise ValueError("The output folder is already populated, please clean up or use another one.")

    print("Started at: ", str(dt.datetime.now()))


    print(str(dt.datetime.now()) + " Reading the GTF file.")
    gtf = read_gtf(args.gtf, subset_kind = args.feature)

    print("Getting the alignment informations")
    aligned_mols, not_passed = get_aligned_molecules(args.bam, gtf, reads = args.reads)

    print(str(not_passed) + " alignments were not associated with a feature.")

    print(str(dt.datetime.now()) + " Counting molecules")
    counted = parse_molecules(aligned_mols)

    print(str(dt.datetime.now()) + " Writing output")
    make_sparse_mtx(counted, args.outdir)

    print(str(dt.datetime.now()) + " Finished")

