import os
import yaml
import sys

if len(sys.argv) == 3:
    cyrius_fname = sys.argv[1]
    if os.path.isfile(cyrius_fname):
        print("Using {} to generate CYP2D6 gt_to_pheno file for this batch.".format(cyrius_fname))
    else:
        print("Cyrius output file does not exist. Exiting.")
        sys.exit(1)

    if os.path.isfile(sys.argv[2]):
        print("Using the following config file: ", sys.argv[2])
        config = yaml.safe_load(open(str(sys.argv[2])))
    else:
        print("Incorrect config file supplied. Exiting.")
        sys.exit(1)
else:
    print("Error: Incorrect number of parameters supplied")
    sys.exit(1)

# Set gt_to_score and gt_to_pheno filenames and paths
gene_anno_dir = config['pgx']['gene_annotations_dir']
out_dir = config['pgx']['cyp2d6_gt_to_pheno']
gt_to_score_fname = os.path.join(gene_anno_dir, 'CYP2D6/CYP2D6_gt_to_score.txt')
if os.path.isfile(gt_to_score_fname):
    print("Using the following gt_to_score file: ", gt_to_score_fname)
else:
    print("GT to score file does not exist. Exiting.")
    sys.exit(1)

gt_to_pheno_fname = os.path.join(out_dir, 'CYP2D6_gt_to_pheno.txt')
print("Path to CYP2D6 gt to pheno file for this batch: {}".format(gt_to_pheno_fname))

# Read gt_to_score file to generate activity score dictionary
with open(gt_to_score_fname, 'r') as f:
    header_fields = next(f).rstrip().split('\t')
    expected_headers = ['Allele', 'Activity Value', 'Source']
    if not(all(expected_header in header_fields for expected_header in expected_headers)):
        print("Error: Incorrect headers found")
        sys.exit(1)

    score_dict = {}

    for line in f:
        fields = line.rstrip('\n').split('\t')
        score_dict[fields[0]] = fields[1]
f.close()


gt_to_pheno_f = open(gt_to_pheno_fname, 'w+')
gt_to_pheno_f.write('\t'.join(["Sample", "Diplotype", "Phenotype", "Diplotype Notes\n"]))

# Read Cyrius output file 
with open(cyrius_fname, 'r') as f:
    header_fields = next(f).rstrip().split('\t')
    expected_headers = ['Sample', 'Genotype', 'Filter']
    if not(all(expected_header in header_fields for expected_header in expected_headers)):
        print("Error: Incorrect headers found")

    gt_pheno_dict = {}

    for line in f:
        fields = line.rstrip('\n').split('\t')
        if fields[2] == 'PASS':
            sample = fields[0]
            geno = fields[1]
            geno_parts = [x.rstrip('+|/') for x in geno.split('*') if x != '']
            score = ""
            for part in geno_parts:
                if 'x' in part:
                    score += ''.join([score_dict[''.join(['*',part[0:part.index('x')]])], '*', part[part.index('x')+1:]]) + '+'
                else :
                    score += score_dict[''.join(['*',part])] + '+'

            k = 'CYP2D6 ' + geno
            pheno = ""

            if 'N/A' in score:
                pheno = "CYP2D6 Indeterminate"
                diplotype_note = "This result indicates that the patient carries at least one uncertain and/or unknown function allele of CYP2D6, and the predicted phenotype for this patient cannot be determined currently for CYP2D6 substrates."
                v = [pheno, diplotype_note+'\n']
                gt_to_pheno_f.write('\t'.join([sample, k, "\t".join(v)]))

            else:
                score = eval(score.rstrip('+'))
                if score == 0:
                    pheno = "CYP2D6 Poor Metabolizer  - Activity score 0"

                elif 0.25 <= score <= 1:
                    pheno = "CYP2D6 Intermediate Metabolizer  - Activity score 0.25 to 1"

                elif 1.25 <= score <= 2.25:
                    pheno = "CYP2D6 Normal Metabolizer  - Activity score 1.25 to 2.25"

                elif score > 2.25:
                    pheno = "CYP2D6 Ultrarapid Metabolizer - Activity score >2.25"

                v = [pheno, '\n']
                gt_to_pheno_f.write('\t'.join([sample, k, "\t".join(v)]))

        elif fields[1] == 'None':
            k = 'CYP2D6 No Call'
            v = ['NA', '\n']
            gt_to_pheno_f.write('\t'.join([sample, k, "\t".join(v)]))

    f.close()

gt_to_pheno_f.close()
