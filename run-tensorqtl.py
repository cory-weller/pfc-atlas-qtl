#!/usr/bin/env python

import argparse
import os
import sys
from pathlib import Path
import logging
import logging.handlers
from tensorqtl import *
import pandas as pd
from datetime import datetime
import statistics
import copy

os.chdir('/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl')

# Set up logging and exit handlers
logger = logging.getLogger('mylogger')
logger.setLevel(logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument('--dry',
                    default=False,
                    action="store_true",
                    help='Validate arguments and file permissions, but do not run')


parser.add_argument('--traw',
                    default=None,
                    type=str,
                    help='Path to Plink transposed raw genotypes file (plain text), where rows correspond to genomic positions. Must contain a single chromosome.')


parser.add_argument('--counts', 
                    default=None,
                    required=True,
                    type=str,
                    help='Path to (pseudobulked) tabular feature counts file.')


parser.add_argument('--counts-index', 
                    default='1',
                    required=False,
                    type=str,
                    help='Column name or number corresponding to Sample ID in the counts file. By default, uses the first column.')


parser.add_argument('--covariates', 
                    default=None,
                    required=True,
                    type=str,
                    help='Path to sample covariates file. All covariates will and samples will be used, when possible.')

parser.add_argument('--covariates-index', 
                    default='1',
                    required=False,
                    type=str,
                    help='Column name or number corresponding to Sample ID in the covariates file. By default, uses the first column.')

parser.add_argument('--interaction', 
                    default=None,
                    required=False,
                    type=str,
                    help='Name of interaction term from covariates file')


parser.add_argument('--outdir', 
                    default='None',
                    required=True,
                    type=str,
                    help='Directory for final output report. Will be created if it does not exist.')


parser.add_argument('--prefix', 
                    default='tensorQTL',
                    type=str,
                    help='String to prefix final report output.')

parser.add_argument('--mode', 
                    default=None,
                    required=True,
                    choices=['rna','atac','RNA','ATAC'],
                    type=str,
                    help='Is count data reflecting gene counts (RNA) or ATAC fragments?')

parser.add_argument('--rna-features', 
                    default='data/GRCh38-2024-A-rna-features-unique.txt',
                    type=str,
                    help='Path to file containing chr, start, end, gene_name matching the counts rna count data.')

parser.add_argument('--window', 
                    default=1000000,
                    type=int,
                    help='')

parser.add_argument('--skiptqtl', 
                    default=False,
                    action='store_true',
                    help='')


# Require --covariates-file
# Optional --covs-include (Default to ALL) OR --covs-exclude (Default to NONE) but NOT both, comma-separated list
# Optional --samples-include (Default to ALL) OR --samples-exclude (Default to None), comma-separated list
# nargs='+' takes 1 or more arguments, nargs='*' takes zero or more.
# Require --counts pointing to file with pseudobulk counts. Autodetect if rows are samples or features?
# Optional --features pointing to a file with list (one per row, no header) 
# Optional --window
# Import test args if none provided

if len(sys.argv) > 1:
    args = parser.parse_args()
    class ExitHandler(logging.StreamHandler):
        def emit(self, record):
            super().emit(record)
            if record.levelno in (logging.ERROR, logging.CRITICAL):
                sys.exit(1)
else:
    import testargs
    logger.info('Using test arguments because none were provided!')
    args = parser.parse_args(testargs.args)
    class ExitHandler(logging.StreamHandler):
        def emit(self, record):
            super().emit(record)
            if record.levelno in (logging.CRITICAL):
                sys.exit(1)


if args.interaction == ['None']:
    args.interaction = None


args.mode = args.mode.upper()



eh = ExitHandler()
logger.addHandler(eh)

# Logging format
logger_format = '%(asctime)s %(module)-15s %(levelname)-15s %(message)s'
date_format = '%Y-%m-%d %H:%M:%S'
formatter = logging.Formatter(logger_format, date_format)
eh.setFormatter(formatter)


if args.interaction is None:
    interaction = 'interaction-None'
else:
    interaction = 'interaction-' + '-'.join(args.interaction)

timestring = datetime.now().strftime('%Y%m%d')  # YYYYMMDD formatted string
args.outdir = str(Path(args.outdir).absolute())     # Get absolute path


def abspath(filepath:str) -> str:
    i = str(Path(filepath).absolute())
    if not os.path.isfile(i):
        logger.warning(f'{i} does not exist or is spelled incorrectly')
    return i


# import torch
# from tensorqtl import pgen, cis, trans, post
# device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas {pd.__version__}")


# Define absolute file paths
trawFile =        abspath(args.traw)
covariatesFile =    abspath(args.covariates)
countsFile =    abspath(args.counts)



# Check files exist
requiredFiles = [trawFile, covariatesFile, countsFile]

if args.mode == 'RNA':
    rnaFeaturesFile = abspath(args.rna_features)
    requiredFiles.append(rnaFeaturesFile)




missingFiles = 0
for file in requiredFiles:
    if not os.path.isfile(file):
        missingFiles += 1
        logger.info(f'{file} does not exist')


if missingFiles > 0:
    logger.error('Quitting due to missing files')




####################################################################################################
# GENOTYPES IMPORT
####################################################################################################
# genotypes_df = pd.read_csv(trawFile, sep='\t', engine='pyarrow')
genotypes_df = pd.read_csv(trawFile, sep='\t', engine='pyarrow')

# By default, plink outputs genotype column headers as {x}_{x}, so if it needs to be fixed,
# uncomment and run the following lines:
# col_rename = {x:x.split('_')[1] for x in genotypes_df.columns[6:]}
# genotypes_df.rename(columns=col_rename, inplace=True)

# Get list of chromosomes in genotype data for subsetting features
chrs_with_genotypes = genotypes_df['CHR'].unique().tolist()
to_drop = ['(C)M','COUNTED','ALT']
genotypes_df.drop(to_drop, axis=1, inplace=True)
genotypes_df.rename(columns={'SNP':'snp','CHR':'chr','POS':'pos'}, inplace=True)
genotypes_df.sort_values(['chr','pos'], inplace=True)

# Set up variants df
variant_df = copy.deepcopy(genotypes_df[['snp','chr','pos']])
variant_df.set_index('snp', inplace=True)
variant_df.rename(columns={'chr':'chrom'}, inplace=True)

genotypes_df.drop(['chr','pos'], axis=1, inplace=True)

genotypes_df.set_index('snp', inplace=True)

if genotypes_df.shape[0] != variant_df.shape[0]:
    logger.error('Unexpected: counts_df and pos_df are different lengths')

# Force all columns to numeric dtype
genotypes_df = genotypes_df.apply(lambda col:pd.to_numeric(col, errors='coerce'))



####################################################################################################
# COUNTS IMPORT
####################################################################################################
# Import as dataframe
# Note: sep=None, engine='python' enables auto-detection of type and delimiter
# counts_df = pd.read_csv(countsFile, sep=',', engine='pyarrow')
counts_df = pd.read_csv(countsFile, sep=',', engine='c')

# Drop 'cohort' column if it exists
# Add more columns to the list if needed
to_drop = ['cohort','cell_type']
counts_df.drop(to_drop, axis=1, errors='ignore', inplace=True)

# Define the column to index on (sample names)
try:
    # Attempt to convert to integer. Only successful if an integer was provided as argument.
    # Otherwise it will produce an ignored error and continue with the string as provided.
    args.counts_index = counts_df.columns[int(args.counts_index) - 1]
except ValueError:
    pass

counts_samples = [x for x in counts_df[args.counts_index]]





# Set up transposed counts table
# rows = 'phenotype_id' genes or ATAC features
# columns = SampleIDs
counts_df = counts_df.transpose()
counts_df.columns  = [x for x in counts_df.loc[args.counts_index]]
counts_df.drop(args.counts_index, inplace=True)

#counts_df = counts_df.convert_dtypes()


# Build annotation 
if args.mode == 'RNA':
    features_df = pd.read_csv(rnaFeaturesFile, sep=None, engine='python')
    # Drop everything except chr1-22
    autosomes = [f'chr{x}' for x in range(1,23)]
    features_df = features_df[features_df['chr'].isin(autosomes)]
    features_df = features_df[features_df['chr'].isin(chrs_with_genotypes)]
    features_df.rename(columns={'gene_name':'phenotype_id'}, inplace=True)
    features_df.set_index('phenotype_id', inplace=True)
    # Merge annotation with counts
    counts_df = pd.merge(counts_df, features_df, left_index=True, right_index=True)
    counts_df.reset_index(inplace=True)
    counts_df.rename(columns={'index':'phenotype_id'}, inplace=True)
elif args.mode == 'ATAC':
    counts_df['chr'] = [x.split(':')[0] for x in counts_df.index.tolist()]
    counts_df = counts_df[counts_df['chr'].isin(chrs_with_genotypes)]
    atacFeatures = counts_df.index.tolist()
    peakPositions = [x.split(':')[1].split('-') for x in atacFeatures]
    counts_df['start'] = [round(statistics.mean([int(x) for x in y])) for y in peakPositions]
    counts_df['end'] = [x+1 for x in list(counts_df['start'])]
    counts_df['phenotype_id'] = [x.replace(':','_').replace('-','_') for x in atacFeatures]
    counts_df.reset_index(inplace=True)
    counts_df.drop('index', axis=1, inplace=True)
    counts_df.rename(columns={'#chr':'chr'}, inplace=True)


# Reorganize columns
counts_df.insert(0, 'end', counts_df.pop('end'))
counts_df.insert(0, 'start', counts_df.pop('start'))
counts_df.insert(0, 'chr', counts_df.pop('chr'))
counts_df.insert(0, 'phenotype_id', counts_df.pop('phenotype_id'))
#counts_df.rename(columns={i:i.lower().replace('#chr','chr') for i in counts_df.columns[:3]}, inplace=True)
counts_df['start'] += 1  # change to 1-based
counts_df['end'] += 1  # change to 1-based


# Sort on chr, start, end
counts_df = counts_df.sort_values(['chr','start','end'])
counts_df = counts_df.set_index('phenotype_id')



pos_df = copy.deepcopy(counts_df[['chr', 'start']])
pos_df.rename(columns={'start':'pos'}, inplace=True)
counts_df.drop(['chr', 'start', 'end'], axis=1, inplace=True)


if pos_df.shape[0] != counts_df.shape[0]:
    logger.error('Unexpected: counts_df and pos_df are different lengths')


# Force to float64 (NOT capital F Float64)
counts_df = counts_df.apply(lambda col:pd.to_numeric(col, errors='coerce'))


counts_df # index=phenotype_id, sample1, sample2, ... sampleN
pos_df  # no Index, chr, pos

genotypes_df # No Index, sample1, sample2, ... sampleN
variant_df  # Index=snp, chr, pos

genotypes_samples = genotypes_df.columns.tolist()



####################################################################################################
# COVARIATES TABLE PREP
####################################################################################################


# Load covariates
covariates_df = pd.read_csv(covariatesFile, sep=None, engine='python')
covariates_samples = covariates_df[args.covariates_index].tolist()
covariates_df.set_index(args.covariates_index, inplace=True)

# Get sample intersect
intersect = [x for x in counts_samples if x in genotypes_samples]
intersect = [x for x in intersect if x in covariates_samples]
intersect = sorted(intersect)


# Subset all tables to intersection of samples
genotype_df_intersect = genotypes_df[intersect]
counts_df_intersect = counts_df[intersect]
covariates_df_intersect = covariates_df[covariates_df.index.isin(intersect)]

# Filter genotypes to 5% MAF
genotypes_df['altsum']=genotypes_df.sum(axis=1)

allele_min = 0.05*2*len(intersect)
allele_max = 0.95*2*len(intersect)
genotypes_df = genotypes_df[genotypes_df.altsum > allele_min]
genotypes_df = genotypes_df[genotypes_df.altsum < allele_max]
genotypes_df.drop('altsum', axis=1, inplace=True)


variant_df = variant_df.loc[genotypes_df.index.tolist()]


# Ensure output directory exists
os.makedirs(args.outdir, exist_ok=True) 


with open(f'{args.outdir}/model_samples.txt', 'w') as outfile:
    outfile.write('\n'.join(intersect)+'\n')


with open(f'{args.outdir}/model_features.txt', 'w') as outfile:
    outfile.write('\n'.join(counts_df.index.tolist())+'\n')


with open(f'{args.outdir}/model_variants.txt', 'w') as outfile:
    outfile.write('\n'.join(genotypes_df.index.tolist())+'\n')


with open(f'{args.outdir}/model_summary.txt', 'w') as outfile:
    outfile.write(f'Samples: {len(intersect)}\n')
    outfile.write(f'Features: {len(counts_df.index.tolist())}\n')
    outfile.write(f'Variants: {len(genotypes_df.index.tolist())}\n')


# # Load modality covariates
# if args.modality_covariates is not None:
#     modality_covariates_df = pd.read_csv(modalityCovariatesFile, sep=None, engine='python')
#     modality_covariates_df = modality_covariates_df.rename(columns={'sample':'FID'})
#     modality_covariates_df.set_index('FID', inplace=True)
#     for term in args.modality_covariates:
#         if term not in modality_covariates_df.columns:
#             logger.error(f"Interaction term {term} is not in provided modality covariates file!")
#     modality_covariates_df = modality_covariates_df[args.modality_covariates]
#     modality_covariates_df.rename(columns=lambda x: f'{args.mode}_{x}', inplace=True)


# # Get interactions
# if args.interaction is not None:
#     for term in args.interaction:
#         if term not in covariates_df.columns:
#             logger.error(f"Interaction term {term} is not in provided covariates file!")
#     interaction_df = covariates_df[args.interaction]


# # Subset to desired covariates
# for term in args.covariates:
#     if term not in covariates_df.columns:
#         logger.error(f"Covariate term {term} is not in provided covariates file!")


# covariates_df = covariates_df[args.covariates]
# covariates_df.rename(columns=lambda x: x.replace('PC', 'geno_PC'), inplace=True)

# # covariates_df = modality_covariates_df.join(covariates_df, how='inner')


# # Get genotypes
# from tensorqtl import genotypeio, cis, trans


# pr = genotypeio.PlinkReader(plink_prefix_path)
# # load genotypes and variants into data frames
# genotype_df = pr.load_genotypes()

# if args.interaction is not None:
#     interaction_df_intersect = interaction_df.loc[interaction_df.index.isin(intersect)]
# else:
#     interaction_df_intersect = None


# Ensure variant_df chromosome ids are formatted as 'chrN' instead of 'N'

if False:
    variant_df['chrom'] = [f'chr{x}' for x in list(variant_df['chrom'])]




# Subset chromosome specified by --chr
# counts_df_intersect = counts_df_intersect.loc[pos_df['chr'] == args.chr]
# pos_df = pos_df.loc[pos_df['chr'] == args.chr]

if args.skiptqtl:
    sys.exit(0)

# Permutations NOT RELLY YET IMPLEMENTED
if False:
    if args.qtlmethod == 'permute':
        print('Doing permutations')
        cis_df = cis.map_cis(genotype_df = genotype_df_intersect,
                                variant_df = variant_df,
                                phenotype_df = counts_df_intersect,
                                phenotype_pos_df = pos_df, 
                                covariates_df = covariates_df_intersect,
                                window=args.window)
        print('Calculating qvalues')
        tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)
        print('Outputting test_out.csv')
        cis_df.to_csv('test_out.csv')


# args.qtlmethod == 'nominal':
cis.map_nominal(genotype_df = genotype_df_intersect, 
                variant_df = variant_df, 
                phenotype_df = counts_df_intersect, 
                phenotype_pos_df = pos_df, 
                prefix = 'cis-nominal',
                output_dir = args.outdir,
                covariates_df = covariates_df_intersect,
                window = args.window,
                run_eigenmt=True,
                write_top=True,
                write_stats=True)

#                 interaction_df = interaction_df_intersect,

print('done!')
sys.exit(0)

