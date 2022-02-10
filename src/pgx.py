import subprocess
import argparse
import yaml
import logging
import sys
from pathlib import Path, PurePath
import pgx_utils
from genotype_variants import _genotype_sample
from discover_variants import _annotate_genotypes


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='PGxTyper: Discover & Annotate Clinically Actionable Pharmacogenomic Interactions')
    parser.add_argument('-b', dest='bam', required=True, help='Path to BAM file')
    parser.add_argument('-s', dest='sample_long_id', required=True, help='Full sample identifier')
    parser.add_argument('-a', dest='accession', required=False, help='Sample accession')
    parser.add_argument('-c', dest='config', required=False, help='confg file')
    parser.add_argument('-o', dest='outdir', required=True, help='Output Directory')
    parser.add_argument('-t', dest='thresholds', required=False, help='thresholds to use for FP detection; "WGS" or "WES"', default='WGS')
    parser.add_argument('--skip-pgx-FDA', dest='do_pgx_FDA', required=False, default=True, action='store_false', help='Skip the FDA version of PGx.')
    parser.add_argument('--skip-pgx-CPIC', dest='do_pgx_CPIC', required=False, default=True, action='store_false', help='Skip the CPIC version of PGx.')
    args = parser.parse_args()

    Path(PurePath(args.outdir)).mkdir(parents=True,exist_ok=True)

    logger = logging.getLogger('PgxTyper')
    logger.setLevel(logging.DEBUG)
    logfile = logging.FileHandler(filename=PurePath(args.outdir).joinpath('pgx.log'),mode='a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logfile.setFormatter(formatter)
    logger.addHandler(logfile)

    #Check arguments
    if args.thresholds not in ['WES', 'WGS']:
        logger.error('False positive thresholds must be either "WES" or "WGS" (default)')
        sys.exit(1)
    assert(Path(args.bam).resolve(strict=True))
    assert(Path(args.outdir).resolve(strict=True))
    assert(Path(args.config).resolve(strict=True))

    if args.accession:
        accession = args.accession
    else:
        accession = args.sample_long_id
    config = yaml.safe_load(open(str(args.config)))

    logger.info("Beginning PGxTyper workflow for {}".format(args.sample_long_id))

    test_codes = {}
    if args.do_pgx_FDA or args.do_pgx_CPIC:
        test_codes['pgx'] = config['pgx']['test_code']
        Path(PurePath(args.outdir)).joinpath(test_codes['pgx'], 'supporting').mkdir(parents=True,exist_ok=True)

    for workflow in test_codes:
        test_code = test_codes[workflow]
        out_path = Path(PurePath(args.outdir)).joinpath(test_codes[workflow])
        out_filename_prefix = args.sample_long_id+'_'+test_code
        outprefix = PurePath(out_path).joinpath('supporting', out_filename_prefix)
        all_bases_vcf = PurePath(out_path).joinpath('supporting', out_filename_prefix+'.allbases.vcf')
        tmpdir = str(PurePath(args.outdir).joinpath(args.sample_long_id,test_code+"tmp"))
        sex = _genotype_sample(args.bam,config,workflow,outprefix,tmpdir)
        genotyping_report = PurePath(out_path).joinpath('supporting', out_filename_prefix+'.genotype.txt')
        genotyping_report_xlsx = PurePath(out_path).joinpath(out_filename_prefix+'.genotype.xlsx')
        _annotate_genotypes(all_bases_vcf,config[workflow]['genotyping']['roi_bed'],sex,genotyping_report,genotyping_report_xlsx,args.thresholds)
        if workflow == 'pgx':
            diplotypes_phenotypes_report = PurePath(out_path).joinpath(out_filename_prefix+'_diplotypes_phenotypes.txt')
            gene_phenotypes, extra_genes_variants = pgx_utils._diplotypes_to_phenotypes(genotyping_report, config['pgx']['gene_annotations_dir'], config['pgx']['genes_without_annotations'], config['pgx']['cyp2d6_gt_to_pheno'], args.accession)
            if args.do_pgx_CPIC:
                phenotypes_CPIC_drug_info, drugs_per_gene = pgx_utils._phenotypes_to_CPIC_drug_info(gene_phenotypes, config['pgx']['gene_annotations_dir'])
                pgx_report_cpic_fname = PurePath(out_path).joinpath(out_filename_prefix+".CPIC_report.xlsx")
                pgx_utils._write_CPIC_output(gene_phenotypes, phenotypes_CPIC_drug_info, drugs_per_gene, config['pgx']['cpic_output_template'], pgx_report_cpic_fname, accession)
            if args.do_pgx_FDA:
                final_report_fda = PurePath(out_path).joinpath(out_filename_prefix+'.FDA_report.xlsx')
                fda_drug_groups = pgx_utils._read_FDA_drug_file(config['pgx']['fda_drugs'])
                fda_drugs_to_output = pgx_utils._gather_FDA_output(gene_phenotypes, fda_drug_groups)
                pgx_utils._write_FDA_output(config['pgx']['fda_output_template'], fda_drugs_to_output, extra_genes_variants, final_report_fda, accession)


