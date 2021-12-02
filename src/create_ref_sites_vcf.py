#Takes a .gvcf file and a snps+indels file, and outputs a vcf file containing non-variant sites plus indels with no annotations and no alt site listed. 
#This file will be merged with the snps+indels file to create an "allbases" vcf file

import sys
import re


def _get_snp_positions(fname):
	positions = set([])
	with opne(fname, 'r') as f:
	for line in f:
		if line[0] != '#': #non-header line
			fields = line.split('\t')
			chr = fields[0]
			pos = fields[1]
			ref = fields[3]
			alt = fields[4]
			if len(ref) == 1 and len(alt) == 1:
				positions.add(chr+':'+pos)
	f.close()
	return(positions)
	
def _create_ref_sites_vcf(gvcf, all_calls_vcf, out_vcf):
	snp_positions = _get_snp_positions(all_calls_vcf)

	with open(out_vcf, 'w') as out:
		with open(gvcf, 'r') as gvcf_file:
			for line in gvcf_file:
				if line[0] == '#':
					out.write(line)
				else:
					fields = line.split('\t')
					chr = fields[0]
					pos = fields[1]
					id = fields[2]
					ref = fields[3]
					alt = fields[4]
					qual = fields[5]
					filter = fields[6]
					info = fields[7]
					format = fields[8]
					sample_data = fields[9]
					dp = sample_data.split(':')[2]
					key = chr+':'+pos
					if key not in snp_positions:
						out.write('\t'.join([chr, pos, id, ref[0], '.', qual,filter,'DP='+dp,format,sample_data]))


