import pandas as pd
import numpy as np
import sys
import re
from subprocess import check_output
from multiprocessing.dummy import Pool as ThreadPool 
import itertools

# create a GTF dict
def read_gtf(filepath, countfeature = "gene_id"):
	gtf_dict = {}
	with open(filepath, "r") as gtffile:
		for entry in gtffile:
			entry = entry.split("\t")

			chrom = entry[0]
			kind = entry[2]
			start = entry[3]
			stop = entry[4]
			strand = entry[6]
			geneID = entry[8]
			geneID = geneID.split(";")[0]
			geneID = geneID.split(" ")[1]
			geneID = re.sub("\"", "", geneID)

			#print(geneID + " " + chrom + " " + kind + " " + start + " " + stop + " " + strand)

			gtf_dict[geneID] = chrom + "," + kind + "," + start + "," + stop + "," + strand
	return(gtf_dict)

def call_samtools(bamfilepath, gtf_dict, geneid, strandness = "sense", alignment_qual = "10"):
	'''call samtools
	the strandness parameter defines whether the gene-sense or antisense reads should be counted
	'''
	chrom, kind, start, stop, strand = gtf_dict[geneid].split(",")

	# Invoking samtools to extract the reads within the specified region
	# A minimal alignment quality of 10 is required
	if strandness == "sense":
		if strand == "+":
			samtools_out = check_output(["samtools", "view", "-q", alignment_qual, "-F 16", bamfilepath, chrom + ":" + start + "-" + stop]).decode("utf-8")
		else:
			samtools_out = check_output(["samtools", "view", "-q", alignment_qual, "-f 16", bamfilepath, chrom + ":" + start + "-" + stop]).decode("utf-8")
	elif strandness == "antisense":
		if strand == "+":
			samtools_out = check_output(["samtools", "view", "-q", alignment_qual, "-f 16", bamfilepath, chrom + ":" + start + "-" + stop]).decode("utf-8")
		else:
			samtools_out = check_output(["samtools", "view", "-q", alignment_qual, "-F 16", bamfilepath, chrom + ":" + start + "-" + stop]).decode("utf-8")

	# Getting the reads within this region
	reads = samtools_out.split("\n") # the -1 removes the last empty entry due to the split at "\n"
	del reads[-1]

	# Processing the retrieved reads
	cell_dict = {}  # used for storing the UMI count per cell-barcode
	for read in reads:
		slots = read.split("\t")
		# Retrieving the barcode sequence

		barcode = [x for x in slots if x.startswith("CB")] # Test whether barcode slot exists
		if len(barcode) > 0:
			barcode = barcode[0]
			barcode = re.sub("CB\:Z\:", "", barcode)

			# Retrieving the UMI sequence
			umi = [x for x in slots if x.startswith("UB")] # Test whether UMI slot exists
			if len(umi) > 0:
				umi = umi[0]
				umi = re.sub("UB\:Z\:", "", umi)

				if barcode in cell_dict.keys():
					if umi not in cell_dict[barcode]:
						cell_dict[barcode].append(umi)
				else:
					cell_dict[barcode] = [umi]


	cellIDs = cell_dict.keys()
	umi_count = 0
	for cellID in cellIDs:
		cell_dict[cellID] = len(cell_dict[cellID])
		#print(geneid + " " + cellID + " " + str(cell_dict[cellID]))
		umi_count += cell_dict[cellID]
	# Getting the number of reads in this region
	testlen = len(reads)

	#print([geneid, testlen, len(cell_dict.keys()), umi_count])

	return([geneid, cell_dict])


if __name__ == "__main__":
	gtf_file_path = sys.argv[1]
	bam_file_path = sys.argv[2]
	strandness = sys.argv[3]
	out_file_path = sys.argv[4]

	gtf_dict = read_gtf(gtf_file_path)
	print("A GTF file with %i features was read" % len(gtf_dict.keys()))

	# Generating a list of geneIDs
	genes = list(gtf_dict.keys())

	# Spawning a pool for multithreading
	with open(out_file_path, "w+") as outfile:
		with ThreadPool(4) as pool:
			for geneID, countdict in pool.starmap(call_samtools, zip(itertools.repeat(bam_file_path), itertools.repeat(gtf_dict), genes, itertools.repeat(strandness))):
				if len(countdict.keys()) > 0:
					for cell in list(countdict.keys()):
						outfile.write(geneID + "\t" + cell + "\t" + str(countdict[cell]) + "\n")

