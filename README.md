# lmPGx v2.1 Pipeline
Performs genotyping, diplotype identification, phenotyping, and reporting

Warning/Disclaimer:
The code and data are provided without warranty and the user accepts all risk and responsibility for testing and validating use of lmPGX in their own environment. lmPGX is not FDA approved, nor has approval been sought from the FDA.

## Setup
### Create a virtual environment
From the src directory:<br/> 
`python3 -m venv env`<br/>
`source env/bin/activate`<br/>
`python3 -m pip install -r requirements.txt`

### Download dependencies and update config file
lmPGX requires the following:
* python3
* GATK4 (tested with GATK4.0.3.0)
* dbsnp VCF file (All_YYYYMMDD.vcf.gz)
* java (tested with 1.8.0
* picard (tested with 2.18.3
* CYP2D6 genotype TSV file (containing the header: Sample\tGenotype\tFilter) generated by Cyrius or a similar tool

Update the paths to these tools/resources in lib/config.yml, as well as the other paths specified by lib/config.yml (reference FASTA file and GATK temp out directory). It's best to update all paths in lib/config.yml to be absolute paths. The path for cyp2d6_gt_to_pheno file will be where the cyp2d6_gt_to_pheno file will be generated by the create_cyp2d6_gt2pheno.py script. This gt_to_pheno file will then be used by the PGx pipeline.

### Example Usage
* activate Python venv
source env/bin/activate

* Generate the genotype to phenotype file for CYP2D6 using the CYP2D6 genotype TSV file
`python3 src/create_cyp2d6_gt2pheno.py CYP2D6_genotype_file.tsv lib/config.yml`

* Run the PGx pipeline
`python3 src/pgx.py -b BAM/CRAM -s SAMPLE_ID  -o OUTDIR -c lib/config.yml`
