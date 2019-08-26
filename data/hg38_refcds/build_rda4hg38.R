library(dndscv)
item=1:22
chr=paste0("chr",item)
path_cds_table = "./hg384dndssv.txt"
path_genome_fasta = "/scratch/cqs/baiy7/Tim_proj/Family_WGS/genome/hg38/Homo_sapiens_assembly38.fasta"
buildref(cdsfile=path_cds_table, genomefile=path_genome_fasta, outfile = "hg38_refcds.rda", onlychrs = chr)
