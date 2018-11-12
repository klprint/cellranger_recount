#!/usr/bin/env python
###
# Created: 26.06.2018
# By: Kevin Leiss
###

import pandas as pd
#import numpy as np
import sys
import time
import re
import itertools
from subprocess import check_output
from multiprocessing import Pool as ThreadPool
import gzip
import parmap
import tqdm

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

                

def call_samtools(bamfilepath, chrom, start, stop, strand, alignment_qual="255", strandness="sense"):
    '''
    Calls samtools for a specific region

    It is important to have samtools in the $PATH variable available!
    bafilepath = path to the indexed bamfile
    chrom = chromosome identifier
    start = start nucleotide
    stop = stop nucleotide
    strand = + or -
    alignment_qual = adjusted for STAR, default = "255" -- needs to be a string
    strandness = returns either reads sense or antisense to the feature

    Returns a list of aligned reads (str)
    '''

    if strandness == "antisense":
        if strand == "+":
            strand = "-"
        else:
            strand = "+"

    if strandness == "sense" or strandness == "antisense":
        if strand == "+":
            aligning_reads = check_output(["samtools", "view", "-q", alignment_qual, "-F 16", bamfilepath, chrom + ":" + start + "-" + stop]).decode("utf-8").split("\n")
        else:
            aligning_reads = check_output(["samtools", "view", "-q", alignment_qual, "-f 16", bamfilepath, chrom + ":" + start + "-" + stop]).decode("utf-8").split("\n")
    else:
        aligning_reads = check_output(["samtools", "view", "-q", alignment_qual, bamfilepath, chrom + ":" + start + "-" + stop]).decode("utf-8").split("\n")


    del aligning_reads[-1]
    return aligning_reads



def get_read_metadata(readstring):
    '''
    Uses a string of a SAM alignment and returns a dictionary of optional flags (last columns of a SAM alignment),
    using the flags as keys and all values as strings.
    '''
    read_info = readstring.split("\t")
    meta_block = read_info[11:len(read_info)]
    meta_block = {re.split(":.:",m)[0]: re.split(":.:",m)[1] for m in meta_block}

    return meta_block

def filter_reads(readlist, cr_altered = False, unique = True):
    '''
    Filters a list of SAM alignments

    Removes all reads which aligned mutliple times and where cellranger gave a "MM" flag.
    Also removes alignments, where no corrected barcode or UMI could be found.
    '''
    out_readslist = list()

    for read in readlist:
        meta_block = get_read_metadata(read)

        if(meta_block["NH"] == "1" and "MM" not in meta_block.keys() and "UB" in meta_block.keys() and "CB" in meta_block.keys()):
            out_readslist.append(read)

    return(out_readslist)

def assign_read_identitiy(read):
    '''
    Parse a SAM read string and return the "CB" and "UB" flag values as tuple.
    '''
    meta_block = get_read_metadata(read)
    return (meta_block["CB"], meta_block["UB"])


def call_samtools_by_gtfdict(gene, bamfilepath, gtfdict, alignment_qual="255", strandness="sense", cr_altered = False, unique = True):
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

    # out = list()

    # with ThreadPool(ncores) as pool:
    #     for out in pool.starmap(call_samtools_by_gtfdict, zip(itertools.repeat(bamfilepath), genes, itertools.repeat(gtf_dict), itertools.repeat(alignment_qual), itertools.repeat(strandness))):
    #         out = out + o
    # i = 1

    # with ThreadPool(ncores) as pool:
    #     out = pool.starmap(call_samtools_by_gtfdict, zip(itertools.repeat(bamfilepath), genes, itertools.repeat(gtf_dict), itertools.repeat(alignment_qual), itertools.repeat(strandness)))


    # with open("testout.txt", "w") as outfile:
    #     for gene in genes:
    #         out = call_samtools_by_gtfdict(bamfilepath, gene, gtf_dict)



    #         for molecule in out:
    #             outfile.write(molecule[0] + "," + molecule[1] + "," + molecule[2] + "\n")


    #         updt(len(genes), i)
    #         i += 1

    out = parmap.map(call_samtools_by_gtfdict, genes, bamfilepath=bamfilepath, gtfdict=gtf_dict, strandness=strandness, pm_pbar=True,  pm_processes=ncores)

    print("Unlisting the outputlists")
    out = list(itertools.chain.from_iterable(out))

    if(umi):
        print("Generating UMI")
        out = list(set(out))
    
    return out


def count_per_bc(gene_bc_umi_tuplelist):
    outdict = dict()

    i = 1
    for entry in gene_bc_umi_tuplelist:
        gene, barcode, umi = entry

        if barcode not in outdict.keys():
            outdict[barcode] = dict()

        if gene not in outdict[barcode].keys():
            outdict[barcode][gene] = 1
        else:
            outdict[barcode][gene] += 1
        updt(len(gene_bc_umi_tuplelist), i)
        i+=1

    return(outdict)
    

def write_sparse_matrix(quant_dict, outdir):
    outdir = re.sub("/", "", outdir)
    barcodes = list(quant_dict.keys())
    genes = set()

    for key, value in quant_dict.items():
        entry_genes = value.keys()
        for gene in entry_genes:
            genes.add(gene)

    genes = list(genes)

    lines = list()
    sum_umi = 0

    for bc, gene_quant in quant_dict.items():
        bc_id = barcodes.index(bc) +1

        for gene, value in gene_quant.items():
            gene_id = genes.index(gene) +1
            sum_umi += value
            value = str(value)

            lines.append(str(gene_id) + " " + str(bc_id) + " " + str(value))



    with open(outdir + "/matrix.mtx", "w") as mtx:
        mtx.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        mtx.write(str(len(genes)) + " " + str(len(barcodes)) + " " +  str(sum_umi) + "\n")
        for line in lines:
            mtx.write(line + "\n")

    with open(outdir + "/genes.tsv", "w") as gns:
        for gene in genes:
            gns.write(gene + "\t" + gene + "\n")

    with open(outdir + "/barcodes.tsv", "w") as bcds:
        for barcode in barcodes:
            bcds.write(barcode + "\n")






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

    parser.add_argument("-s", "--strandness",
        metavar="SENSE|ANTISENSE|IGNORE",
        type=str,
        default="sense",
        help="Define whether only reads in FEATURE sense / antisense direction should be counted, or strand-unspecific")


    parser.add_argument("--reads", 
        action = "store_true",
        help = "If set, the returned matrix will be readcounts instead of UMI counts")


    parser.add_argument("--ncores",
        metavar = "INT",
        type=int,
        default=4,
        help="Number of cores to alocate")


    args = parser.parse_args()

    #print(args)


    import datetime as dt
    import random
    import os

    logfile = "log_" + re.sub(" |:|\\.", "_", str(dt.datetime.now()))[0:19] + ".txt"
    sys.stdout = open(logfile, 'w',  buffering = 1)

    print("StdOut will be redirected into file: " + logfile)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    else:
        if os.path.exists(args.outdir + "/matrix.mtx") or os.path.exists(args.outdir + "/genes.tsv") or os.path.exists(args.outdir + "/barcodes.tsv"):
            raise ValueError("The output folder is already populated, please clean up or use another one.")

    print(str(dt.datetime.now()) + " starting")


    print(str(dt.datetime.now()) + " Reading the GTF file.")
    gtf = read_gtf(args.gtf, subset_kind = args.feature)

    print()
    print(10 * "-")
    print()

    print(str(dt.datetime.now()) + " Reading alignments for " + str(len(gtf.keys())) + " genes")
    molecule_info = check_reads_per_feature(args.bam, gtf, ncores = args.ncores, umi = not args.reads, strandness = args.strandness.tolower())
    print(str(dt.datetime.now()) + " Done")

    print()
    print(10 * "-")
    print()

    print(str(dt.datetime.now()) + " Parsing the data")
    counts_dict = count_per_bc(molecule_info)
    print(str(dt.datetime.now()) + " Done")

    print()
    print(10 * "-")
    print()

    print(str(dt.datetime.now()) + " Writing the output")
    write_sparse_matrix(counts_dict, args.outdir)
    print(str(dt.datetime.now()) + " Done")

    print()
    print(10 * "-")
    print()

    print(str(dt.datetime.now()) + " Finished")

