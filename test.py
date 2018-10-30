import cellranger_recount as cr

gtf_file = cr.read_gtf("test.gtf")
print(gtf_file)

reads = cr.dcall_samtools("mnt/outs/possorted_genome_bam.bam", "1", "6229959", "6230073", "+")
print(reads)
print("----")

p_read = [cr.parse_read(x) for x in reads]

print(p_read)