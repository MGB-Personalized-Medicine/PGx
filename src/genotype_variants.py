import subprocess
import logging
from pathlib import Path, PurePath

genotype_variants_logger = logging.getLogger('PgxTyper.GenotypeVariants')


def _get_snp_positions(fname):
	positions = set([])
	with open(fname, 'r') as f:
		for line in f:
			if line[0] != '#': #non-header line
				fields = line.split('\t')
				chr = fields[0]
				pos = fields[1]
				ref = fields[3]
				alt = fields[4]
				if len(ref) == 1 and len(alt) == 1:
					positions.add(chr+':'+pos)
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

def _genotype_sample(bam,config,workflow,outprefix,tmpdir):
	gvcf = str(outprefix) + '.g.vcf'
	all_calls_vcf = str(outprefix) + '.allcalls.vcf'
	ref_positions_vcf = str(outprefix) + '.ref_positions.vcf'
	all_bases_vcf = str(outprefix) + '.allbases.vcf'

	#1) HaplotypeCaller
	haplotype_caller_command = [
		config[workflow]['genotyping']['gatk'], '--java-options', config[workflow]['genotyping']['java_options'],'HaplotypeCaller',
		'--input', bam,
		'--output', gvcf,
		'--reference', config[workflow]['genotyping']['reference_fasta'],
		'--intervals', config[workflow]['genotyping']['roi_bed'],
		'--dbsnp', config[workflow]['genotyping']['dbsnp'],
		'-mbq', '10',
		'-stand-call-conf', '30',
		'-A', 'QualByDepth',
		'-A', 'FisherStrand',
		'-A', 'RMSMappingQuality',
		'-A', 'MappingQualityZero',
		'-A', 'StrandOddsRatio',
		'-A', 'DepthPerAlleleBySample',
		'--read-filter', 'MappingQualityReadFilter',
		'--minimum-mapping-quality', '17', 
		'--read-filter', 'MappingQualityNotZeroReadFilter',
		'-ERC', 'BP_RESOLUTION']
	genotype_variants_logger.info("Running HaplotypeCaller")
	run_haplotype_caller = subprocess.Popen(haplotype_caller_command,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = run_haplotype_caller.communicate()
	genotype_variants_logger.info(stderr)

	#2) GenotypeGVCFs
	genotype_gvcfs_command = [
		config[workflow]['genotyping']['gatk'], '--java-options', config[workflow]['genotyping']['java_options'], 'GenotypeGVCFs',
		'--variant', str(gvcf),
		'--output', str(all_calls_vcf),
		'--reference', config[workflow]['genotyping']['reference_fasta'],
		'--intervals', config[workflow]['genotyping']['roi_bed'],
		'--dbsnp', config[workflow]['genotyping']['dbsnp'],
		'-A', 'QualByDepth',
		'-A', 'FisherStrand',
		'-A', 'RMSMappingQuality',
		'-A', 'MappingQualityZero',
		'-A', 'StrandOddsRatio',
		'-A', 'DepthPerAlleleBySample',
		'--read-filter', 'MappingQualityReadFilter',
		'--minimum-mapping-quality', '17', 
		'--read-filter', 'MappingQualityNotZeroReadFilter']
	genotype_variants_logger.info("Genotyping GVCF")	
	genotype_gvcf = subprocess.Popen(genotype_gvcfs_command,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = genotype_gvcf.communicate()
	genotype_variants_logger.info(stderr)

	#3) Create REF Sites VCF
	genotype_variants_logger.info("Creating REF Sites VCF")
	_create_ref_sites_vcf(gvcf, all_calls_vcf, ref_positions_vcf)
	
	#4) SortVcf (to combine into allbases.vcf)
	sort_vcf_command = [config[workflow]['genotyping']['java']]
	sort_vcf_command.extend(config[workflow]['genotyping']['java_options'].split(' '))
	sort_vcf_command.extend(['-jar', config[workflow]['genotyping']['picard'], 'SortVcf',
		'I='+ref_positions_vcf,
		'I='+all_calls_vcf,
		'O='+all_bases_vcf])

	genotype_variants_logger.info("Creating All Bases VCF")
	sort_vcf = subprocess.Popen(sort_vcf_command,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = sort_vcf.communicate()
	genotype_variants_logger.info(stderr)

	assert(Path(all_bases_vcf).resolve(strict=True)) 

	sex = 'F'
	return sex
	
