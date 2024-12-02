#!/usr/bin/env python

import argparse
import os
import sys
from pathlib import Path
import logging
import logging.handlers
import tensorqtl
import pandas as pd
from datetime import datetime
import statistics

os.chdir('/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl')

# Set up logging and exit handlers
logger = logging.getLogger('mylogger')
logger.setLevel(logging.DEBUG)



parser = argparse.ArgumentParser()
parser.add_argument('--dry',
                    default=False,
                    action="store_true",
                    help='Validate arguments and file permissions, but do not run DIA-NN')
# parser.add_argument('--plinkfile', 
#                     default=None,
#                     required=True,
#                     type=str,
#                     help='Path to PLINK2 file set {bed,bim,fam}')
parser.add_argument('--celltype', 
                    default=None,
                    required=True,
                    type=str,
                    choices=['Astro', 'ExN', 'InN', 'MG', 'Oligo', 'OPC', 'VC'],
                    help='Celltype')
parser.add_argument('--chr', 
                    default=None,
                    required=True,
                    type=str,
                    choices=[f'chr{x}' for x in list(range(1,23))],
                    help='Celltype')
parser.add_argument('--mode', 
                    default=None,
                    required=True,
                    type=str,
                    choices=['rna', 'atac'],
                    help='Mode')
parser.add_argument('--cohort', 
                    default=None,
                    required=True,
                    type=str,
                    choices=['HBCC', 'NABEC'],
                    help='Cohort')
parser.add_argument('--covariates', 
                    default=None,
                    required=False,
                    type=str,
                    nargs='*',
                    help='List of model covariates to be included from covariates file, e.g. --covariates PC1 PC2 Sex Age')
parser.add_argument('--interaction', 
                    default=None,
                    required=False,
                    type=str,
                    nargs='*',
                    help='Name of interaction term(s) from covariates file')
parser.add_argument('--outdir', 
                    default='None',
                    required=True,
                    type=str,
                    help='File name for final output report. Other output will share the same stem.')

# parser.add_argument('--prefix', 
#                     default='TQTL',
#                     type=str,
#                     help='File name for final output report. Other output will share the same stem.')

parser.add_argument('--qtlmethod', 
                    default=None,
                    type=str,
                    choices=['nominal', 'interaction'],
                    help='File name for final output report. Other output will share the same stem.')
parser.add_argument('--window', 
                    default=1000000,
                    type=int,
                    help='')
parser.add_argument('--pseudobulkmethod', 
                    default=None,
                    type=str,
                    choices=['sum','mean'],
                    help='Pseudobulk / aggregation method')

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

if args.cohort == 'NABEC':
    args.plinkfile = '/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/genotypes/NABEC_polarized.bed'
elif args.cohort == 'HBCC':
    args.plinkfile = '/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl/genotypes/HBCC_polarized_nooutliers.bed'

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
methodsbatch = f'{args.pseudobulkmethod}-{args.qtlmethod}-{timestring}' # level with variable analysis types
countsbatch = f'{args.cohort}-{args.celltype}-{args.mode}-{interaction}'  # Level of celltype counts per mode
args.outdir = str(Path(args.outdir).absolute())     # Get absolute path
args.outdir = f'{args.outdir}/{methodsbatch}/{countsbatch}' # define final output dir for batch

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



# If user provided .bed extension for plink file, remove .bed extension
if args.plinkfile.endswith('.bed'):
    args.plinkfile = args.plinkfile.removesuffix('.bed')

# Define absolute file paths
bedFile =           abspath(f'{args.plinkfile}.bed')
bimFile =           abspath(f'{args.plinkfile}.bim')
famFile =           abspath(f'{args.plinkfile}.fam')
covariatesFile =    abspath(f'data/{args.cohort}-covariates.tsv')
rnaFeatureFile =    abspath('data/rna/RNA-features.tsv')

#countsFile = f'/data/CARD_singlecell/cortex_subtype/output/{args.mode}/{args.celltype}/pseudobulk/{args.cohort.lower()}_cpm_{args.pseudobulkmethod}_log_pseudobulk_model_counts.csv'
countsFile = f'/data/CARD_singlecell/cortex_subtype/output/{args.mode}/{args.celltype}/pseudobulk/{args.cohort.lower()}_cpm_{args.pseudobulkmethod}_pseudobulk_model_counts.csv'

# Check files exist
requiredFiles = [bedFile, bimFile, famFile, countsFile, covariatesFile, rnaFeatureFile]

missingFiles = 0
for file in requiredFiles:
    if not os.path.isfile(file):
        missingFiles += 1
        logger.info(f'{file} does not exist')


if missingFiles > 0:
    logger.error('Quitting due to missing files')


# define path to suffix-removed plink dataset
plink_prefix_path = bedFile.removesuffix('.bed')


# Load and format counts (phenotypes)
# Note: sep=None, engine='python' enables auto-detection of type and delimiter
counts_df = pd.read_csv(countsFile, sep=None, engine='python')  
counts_df.rename(columns={'Unnamed: 0' : 'phenotype_id'}, inplace=True)




# Fix sample ID naming
# NABEC samples remove trailling "-ARC" string
# HBCC samples prepend "HBCC_" and remove trailing "-ARC" string
if args.mode == 'rna':
    firstField = 'phenotype_id'
elif args.mode == 'atac':
    firstField = 'peaks'


sampleids = list(counts_df.columns[1:])
if args.cohort == 'NABEC':
    tokeep = [firstField] + [x for x in sampleids if x.startswith(('UMARY','KEN','SH'))]
    newnames = [firstField] + [f"{x.removesuffix('-ARC')}" for x in tokeep[1:]]
elif args.cohort == 'HBCC':
    tokeep = [firstField] + [x for x in sampleids if not  x.startswith(('U','K','S'))]
    newnames = [firstField] + [f"HBCC_{x.removesuffix('-ARC')}" for x in tokeep[1:]]


counts_df = counts_df[tokeep]
counts_df = counts_df.rename(columns=dict(zip(tokeep, newnames)))
sampleids = list(counts_df.columns[1:])


# Build annotation 
if args.mode == 'atac':
    atacFeatures = list(counts_df['peaks'])
    counts_df['#chr'] = [x.split(':')[0] for x in atacFeatures]
    peakPositions = [x.split(':')[1:3] for x in atacFeatures]
    counts_df['start'] = [round(statistics.mean([int(x) for x in y])) for y in peakPositions]
    counts_df['end'] = [x+1 for x in list(counts_df['start'])]
    counts_df['phenotype_id'] = [x.replace(':','_') for x in list(counts_df['peaks'])]
    counts_df.drop('peaks', axis=1, inplace=True)
    # get chr, start, end from peak column
elif args.mode == 'rna':
    rnaFeatures = pd.read_csv(rnaFeatureFile, sep=None, engine='python')
    counts_df = pd.merge(counts_df, rnaFeatures, on='phenotype_id')


counts_df.insert(0, 'end', counts_df.pop('end'))
counts_df.insert(0, 'start', counts_df.pop('start'))
counts_df.insert(0, '#chr', counts_df.pop('#chr'))
counts_df.insert(0, 'phenotype_id', counts_df.pop('phenotype_id'))
counts_df.rename(columns={i:i.lower().replace('#chr','chr') for i in counts_df.columns[:3]}, inplace=True)
counts_df['start'] += 1  # change to 1-based

# Subset autosomes
allowed_chrs = [f'chr{x}' for x in range(1,23)] # chrs 1-22
counts_df = counts_df.loc[counts_df.chr.isin(allowed_chrs)]

counts_df = counts_df.sort_values(['chr','start','end'])
counts_df = counts_df.set_index('phenotype_id')

pos_df = counts_df[['chr', 'start', 'end']]
counts_df.drop(['chr', 'start', 'end'], axis=1, inplace=True)


# make sure BED file is properly sorted
if not pos_df.equals(pos_df.groupby('chr', sort=False, group_keys=False).apply(lambda x: x.sort_values(['start', 'end']))):
    logger.error("Positions in BED file must be sorted.")

# make sure positions are as expected
if (pos_df['start'] == pos_df['end']).all():
    pos_df = pos_df[['chr', 'end']].rename(columns={'end':'pos'})
else:
    logger.error('Unexpected: start and end positions do not match')

# Load covariates
covariates_df = pd.read_csv(covariatesFile, sep=None, engine='python')
covariates_df.set_index('FID', inplace=True)


# Get interactions
if args.interaction is not None:
    for term in args.interaction:
        if term not in covariates_df.columns:
            logger.error(f"Interaction term {term} is not in provided covariates file!")
    interaction_df = covariates_df[args.interaction]


# Subset to desired covariates
for term in args.covariates:
    if term not in covariates_df.columns:
        logger.error(f"Covariate term {term} is not in provided covariates file!")

covariates_df = covariates_df[args.covariates]


# Get genotypes
from tensorqtl import genotypeio, cis, trans

pr = genotypeio.PlinkReader(plink_prefix_path)
# load genotypes and variants into data frames
genotype_df = pr.load_genotypes()


variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# Get list of samples that intersect all data modes
intersect = [x for x in list(genotype_df.columns) if x in list(covariates_df.index)]
intersect = [x for x in intersect if x in list(counts_df.columns)]
if args.interaction is not None:
    intersect = [x for x in intersect if x in list(interaction_df.index)]

intersect = sorted(intersect)

# Not needed?
# counts_df = counts_df.rename(columns=dict(zip(list(counts_df.columns), [x.replace('-','_') for x in list(counts_df.columns)])))

genotype_df_intersect = genotype_df[intersect]
counts_df_intersect = counts_df[intersect]
covariates_df_intersect = covariates_df.loc[covariates_df.index.isin(intersect)]

if args.interaction is not None:
    interaction_df_intersect = interaction_df.loc[interaction_df.index.isin(intersect)]
else:
    interaction_df_intersect = None


# Ensure variant_df chromosome ids are formatted as 'chrN' instead of 'N'
variant_df['chrom'] = [f'chr{x}' for x in list(variant_df['chrom'])]

# Ensure output directory exists
os.makedirs(args.outdir, exist_ok=True) 

# Subset chromosome specified by --chr
counts_df_intersect = counts_df_intersect.loc[pos_df['chr'] == args.chr]
pos_df = pos_df.loc[pos_df['chr'] == args.chr]

# Permutations NOT RELLY YET IMPLEMENTED
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
elif args.qtlmethod == 'nominal':
    cis.map_nominal(genotype_df = genotype_df_intersect, 
                    variant_df = variant_df, 
                    phenotype_df = counts_df_intersect, 
                    phenotype_pos_df = pos_df, 
                    prefix = 'cis-nominal',
                    output_dir = args.outdir,
                    covariates_df = covariates_df_intersect,
                    interaction_df = interaction_df_intersect,
                    window = args.window,
                    run_eigenmt=True,
                    write_top=True,
                    write_stats=True)


print('done!')
sys.exit(0)

