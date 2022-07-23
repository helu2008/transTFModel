import sys
import pandas as pd
import numpy as np
from sklearn.linear_model import LassoCV
from sklearn.metrics import r2_score
from collections import defaultdict
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import fdrcorrection
from functools import reduce
from random import sample
from sklearn.linear_model import LinearRegression
from argparse import ArgumentParser
import timeit
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning
import allel

simplefilter("ignore", category=ConvergenceWarning)
# general params
CUTOFF_PVALUE = 0.05
NUM_SHUFFLES = 100
NUM_ROBUST_ITER = 10
FRACTION_ROBUST_SAMPLES = 0.9

# read the names of input and output files
parser = ArgumentParser(description='Get parameters')
parser.add_argument('-e', dest='eFile', type=str,
                    help='Normalized gene expression values')
parser.add_argument('-p', dest='pFile', type=str,
                    help='Result from Predixcan')
parser.add_argument('-g', dest='gFile', type=str,
                    help='GTEx filtered genotype file')
parser.add_argument('-t', dest='TF2pos', type=str,
                    help='gene_name_to_rsids_before_pph file')
parser.add_argument('-c', dest='covFile', type=str,
                    help='Tissue covariance file')
parser.add_argument('-o', dest='oFile', type=str,
                    help='Output file')

if len(sys.argv) < 13:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

# Gene expression value
gene_expression_path = args.eFile
start = timeit.default_timer()
# remove covariates
cov = pd.read_table(args.covFile, index_col=0)
covariates = cov.iloc[list(range(3))+list(range(5, 20))+[-2, -1], :]


def remove_covariate_component(expression, covariates):
    expression['gene_id'] = expression['gene_id'].apply(lambda x: x.split('.')[0])
    expression = expression.iloc[:, 3:].set_index('gene_id')
    common_subject = list(set(expression.columns).intersection(covariates.columns))
    if common_subject:
        new_expression = expression[common_subject]
        new_covariates = covariates[common_subject]
        y = new_expression.values.T
        X = new_covariates.values.T
        lm = LinearRegression()
        lm.fit(X, y)
        result = y - lm.predict(X)
        return pd.DataFrame(result.T, index=new_expression.index, columns=new_expression.columns)

    
expression = pd.read_table(gene_expression_path, chunksize=100)
chunk_residuals = []
for i, chunk in enumerate(expression):
    new = remove_covariate_component(chunk, covariates)
    chunk_residuals.append(new)
residuals = pd.concat(chunk_residuals).reset_index()

# Cov-residuals gene expression and the result of predixcan
prediction_path = args.pFile
prediction = pd.read_table(prediction_path)
cols = [x.split('.')[0] for x in list(prediction.columns)]
prediction.columns = cols
predixcan = prediction.set_index('IID').T
residuals = residuals.set_index('gene_id')
common_col = list(set((residuals.columns).intersection(predixcan.columns)))
residuals1 = residuals[common_col]
predixcan1 = predixcan[common_col]
common_gene = list(set((residuals1.index).intersection(predixcan1.index)))
residuals2 = residuals1[residuals1.index.isin(common_gene)]
predixcan2 = predixcan1[predixcan1.index.isin(common_gene)]
predixcan2 = predixcan2.loc[residuals2.index, :]

normalized = residuals2
predixcan = predixcan2
diff = normalized - predixcan

full_expression = pd.read_csv(gene_expression_path, sep='\t')
full_expression.gene_id = full_expression.gene_id.apply(lambda x: x.split('.')[0])
full_expression = full_expression.iloc[:, 3:]
full_expression = full_expression.set_index('gene_id')
expression = full_expression[diff.columns]

diff = diff.reset_index()
end = timeit.default_timer()
print('Finished removing covariates. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

# Mapping between gene codes and gene name
ENSG_to_name_file = pd.read_csv('ENSG_to_name.txt')
ENSG_to_name = dict(zip(ENSG_to_name_file['Gene stable ID'], ENSG_to_name_file['HGNC symbol']))
gene_name_to_ENSG = dict(zip(ENSG_to_name_file['HGNC symbol'], ENSG_to_name_file['Gene stable ID']))

# Mapping from gene to its TFs
tf2gene = pd.read_table('TF2Gene.txt')
gene2tf = dict(tf2gene.groupby('Gene').apply(lambda x: list(x.TF)))

# Extract genotype data

TF2pos = defaultdict(set)
TF2rsid = defaultdict(set)
with open(args.TF2pos, 'r') as f:
    line = f.readline()  # header
    for i, line in enumerate(f):
        (TF,rsid, position) = line.strip().split('\t')
        TF2pos[TF].add(position)
        TF2rsid[TF].add(rsid)

ENSG_to_tf_variant_pos = defaultdict(set)
ENSG_to_tf_rsid = defaultdict(set)
for gene in gene2tf:
    if gene in gene_name_to_ENSG:
        ENSG = gene_name_to_ENSG[gene]
        for TF in gene2tf[gene]:
            for pos in TF2pos[TF]:
                ENSG_to_tf_variant_pos[ENSG].add(int(pos))
            for rsid in TF2rsid[TF]:
                ENSG_to_tf_rsid[ENSG].add(rsid)

end = timeit.default_timer()
print('Finished mapping genes to variants. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

# Calculate the R2 scores and select genes with R2 scores equal to or large than 0
r2_scores = {}
for i, row in residuals2.iterrows():
    y = row.values.reshape(-1, 1)
    y_p = predixcan2.loc[i, :].values.reshape(-1, 1)
    lm = LinearRegression()
    lm.fit(y_p, y)
    y_new = lm.predict(y_p)
    r2 = 1-((y_new-y)**2).sum()/((y-y.mean())**2).sum()
    r2_scores[i] = r2
r2_score_csv = pd.DataFrame(r2_scores.items(), columns=['gene', 'r2_score'])
r2_scores = dict(zip(r2_score_csv['gene'], r2_score_csv['r2_score']))
genes = list(dict(filter(lambda x: x[1] >= 0, r2_scores.items())).keys())
genes_with_2_100 = [x for x in genes if len(ENSG_to_tf_variant_pos.get(x,set())) <= 100 and len(ENSG_to_tf_variant_pos.get(x,set())) >= 2]
genes = genes_with_2_100

end = timeit.default_timer()
print('Finished calculating R2. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

# Read VCF and transform to doses
callset = allel.read_vcf(args.gFile)
samples=callset['samples']
pos=callset['variants/POS']
chromosomes=callset['variants/CHROM']
gt = allel.GenotypeArray(callset['calldata/GT'])
all_doses = gt.to_n_alt()

# keep only selected samples
col_indexes = {}
index = 0
for s in samples:
    col_indexes[s] = index
    index = index + 1

selected_sample_index = []
for s in diff.columns:
    if s in col_indexes:
        selected_sample_index.append(col_indexes[s])

doses = all_doses[:, selected_sample_index]
# Build model
nssnp_model = {}
r2score_dict = {}
for gene in genes:
    if gene in ENSG_to_tf_variant_pos:
        X = doses[list(ENSG_to_tf_variant_pos[gene]), :].transpose()
        X = np.where(np.isnan(X), 0, X)
        y = diff[diff.gene_id == gene].values[0, 1:].astype('d')
        lassocv = LassoCV(cv=10, random_state=0)
        lassocv.fit(X,y)
        y_pred = lassocv.predict(X)
        r2 = r2_score(y, y_pred)
        r2score_dict[gene] = r2
        nssnp_model[gene] = lassocv
 
end = timeit.default_timer()
print('Finished building models. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

def shuffle_y_build_model(gene):
    shuffle_list = []
    if gene in ENSG_to_tf_variant_pos and len(ENSG_to_tf_variant_pos[gene]) > 0:
        X = doses[list(ENSG_to_tf_variant_pos[gene]), :].transpose()
        X = np.where(np.isnan(X), 0, X)
        y = diff[diff.gene_id == gene].values[0, 1:].astype('d')
        for i in range(NUM_SHUFFLES):
            y_shuffle = np.random.permutation(y)
            lassocv = LassoCV(cv=10, random_state=0)
            lassocv.fit(X, y_shuffle)
            y_pred = lassocv.predict(X)
            r2 = r2_score(y_shuffle, y_pred)
            shuffle_list.append(r2)
        return shuffle_list

# permute y
p_values_bckgrnd_1 = defaultdict(float)
for gene in genes:
    if gene in ENSG_to_tf_variant_pos:
        shuffle_list = shuffle_y_build_model(gene)
        _, ztest_pvalue = ztest(x1=np.array([r2score_dict[gene]]), x2=np.array(shuffle_list), alternative='larger')
        if np.isnan(ztest_pvalue):
            ztest_pvalue = 1
        p_values_bckgrnd_1[gene] = ztest_pvalue

_, p_value_corrected = fdrcorrection(np.asarray(list(p_values_bckgrnd_1.values())))
p_values_bckgrnd_1_corrected = dict(zip(p_values_bckgrnd_1.keys(), p_value_corrected))
genes_bgrnd_1 = {k: v for k, v in sorted(p_values_bckgrnd_1_corrected.items(), key=lambda x: x[0]) if v < CUTOFF_PVALUE}

end = timeit.default_timer()
print('Finished shuffling y. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

var_set = reduce(lambda a, b: a.union(b), ENSG_to_tf_variant_pos.values())
def generate_shuffled_nssnp():
    shuffled_var = {}
    for gene, v in ENSG_to_tf_variant_pos.items():
        if len(ENSG_to_tf_variant_pos[gene]) > 0:
            shuffled_var[gene] = sample(list(var_set.difference(v)), len(v))
    return shuffled_var

def shuffle_nssnp_build_model(gene):
    shuffle_list = []
    if gene in ENSG_to_tf_variant_pos:
        for j in range(NUM_SHUFFLES):
            shuffled_var = generate_shuffled_nssnp()
            X = doses[list(shuffled_var[gene]), :].transpose()
            X = np.where(np.isnan(X), 0, X)
            y = diff[diff.gene_id == gene].values[0, 1:].astype('d')
            lassocv = LassoCV(cv=10, random_state=0)
            lassocv.fit(X,y)
            y_pred = lassocv.predict(X)
            r2 = r2_score(y, y_pred)
            shuffle_list.append(r2)
        return shuffle_list
    
# permute nssnp
p_values_bckgrnd_2 = defaultdict(float)
for gene in genes:
    if gene in genes_bgrnd_1:
        if gene in ENSG_to_tf_variant_pos and len(ENSG_to_tf_variant_pos[gene]) > 0:
            shuffle_list = shuffle_nssnp_build_model(gene)
            _, ztest_pvalue = ztest(x1=np.array([r2score_dict[gene]]), x2=np.array(shuffle_list), alternative='larger')
            if np.isnan(ztest_pvalue):
                ztest_pvalue = 1
            p_values_bckgrnd_2[gene] = ztest_pvalue

_, p_value_corrected2 = fdrcorrection(np.asarray(list(p_values_bckgrnd_2.values())))
p_value_bckgrnd_2_corrected = dict(zip(p_values_bckgrnd_2.keys(), p_value_corrected2))
genes_after_shuffle_tf = {k: v for k, v in sorted(p_value_bckgrnd_2_corrected.items(), key=lambda x: x[0]) if v < CUTOFF_PVALUE}

end = timeit.default_timer()
print('Finished shuffling TFs. ' + str(int(end-start)))
start = timeit.default_timer()

significant_gene_list = []
for k in genes_after_shuffle_tf.keys():
    significant_gene_list.append(ENSG_to_name[k])

genes = list(genes_after_shuffle_tf.keys())
snp_weight = {}
with open(args.oFile, 'w') as f:
    f.write('Gene\tR2 score of predixcan\tR2 score of our model\n')
    for ensg in genes_after_shuffle_tf.keys():
        f.write(f'{ENSG_to_name[ensg]}\t{r2_scores[ensg]}\t{r2score_dict[ensg]}\n')
    # write the TF SNPs in the model
    f.write('----------------------------\nTF and Weights: \n')
    f.write('Gene\tvariants\n')
    significant_gene_tf_nssnp = {}
    for ensg in nssnp_model.keys():
        model = nssnp_model[ensg]
        gene = ENSG_to_name[ensg]
        variants = list(ENSG_to_tf_rsid[ensg])
        rsid_weights_dict = dict(zip(variants, model.coef_))
        f.write('%s\t%s\n' % (gene, rsid_weights_dict))

# Robustness check
genes_all_samples = genes
random90hitgenes = []
for i in range(NUM_ROBUST_ITER):
    print("Robustness check #" + i)
    genes = genes_all_samples
    diff = normalized - predixcan
    diff = diff.reset_index()
    all_samples = diff.columns[1:]
    sample_length = len(all_samples)
    random_length = int(FRACTION_ROBUST_SAMPLES * sample_length)
    random_col = sample(list(all_samples), random_length)
    diff = diff.loc[:, ['gene_id'] + list(random_col)]
    selected_sample_index = []
    for s in diff.columns:
        if s in col_indexes:
            selected_sample_index.append(col_indexes[s])

    doses = all_doses[:, selected_sample_index]

    nssnp_model = {}
    r2score_dict = {}
    # modeling the gene expression for the sub-sampled one
    for gene in genes:
        if gene in ENSG_to_tf_variant_pos:
            X = doses[list(ENSG_to_tf_variant_pos[gene]), :].transpose()
            X = np.where(np.isnan(X), 0, X)
            y = diff[diff.gene_id==gene].values[0, 1:].astype('d')
            lassocv = LassoCV(cv=10, random_state=0)
            lassocv.fit(X,y)
            y_pred = lassocv.predict(X)
            r2 = r2_score(y, y_pred)
            r2score_dict[gene] = r2
            nssnp_model[gene] = lassocv

    # permute y
    p_value = defaultdict(float)
    for gene in genes:
        if gene in ENSG_to_tf_variant_pos:
            shuffle_list = shuffle_y_build_model(gene)
            _, ztest_pvalue = ztest(x1=np.array([r2score_dict[gene]]), x2=np.array(shuffle_list), alternative='larger')
            if np.isnan(ztest_pvalue):
                ztest_pvalue = 1
            p_value[gene] = ztest_pvalue
    
    _, p_value_corrected = fdrcorrection(np.asarray(list(p_value.values())))
    robust_p_value_bckgrnd_1 = dict(zip(p_value.keys(), p_value_corrected))
    robust_genes_after_shuffle = {k: v for k, v in sorted(robust_p_value_bckgrnd_1.items(), key=lambda x: x[0]) if v < CUTOFF_PVALUE}
    
    # permute nssnp
    p_value2 = defaultdict(float)
    for gene in genes:
        if gene in robust_genes_after_shuffle:
            if gene in ENSG_to_tf_variant_pos and len(ENSG_to_tf_variant_pos[gene]) > 0:
                shuffle_list = shuffle_nssnp_build_model(gene)
                _, ztest_pvalue = ztest(x1=np.array([r2score_dict[gene]]), x2=np.array(shuffle_list), alternative='larger')
                if np.isnan(ztest_pvalue):
                    ztest_pvalue = 1
                p_value2[gene] = ztest_pvalue
    
    _, p_value_corrected2 = fdrcorrection(np.asarray(list(p_value2.values())))
    robust_p_value_bckgrnd_2 = dict(zip(p_value2.keys(), p_value_corrected2))
    robust_genes_after_tf_shuffle = {k: v for k, v in sorted(robust_p_value_bckgrnd_2.items(), key=lambda x: x[0]) if v < CUTOFF_PVALUE}

    significant_genes = [ENSG_to_name[x] for x in list(robust_genes_after_tf_shuffle)]
    random90hitgenes.append(significant_genes)

end = timeit.default_timer()
print('Finished robustness check. ' + str(int(end-start)))
start = timeit.default_timer()

with open(args.oFile, 'a+') as f:
    f.write("--------------\nRobustness check\n")
    for value in random90hitgenes:
        f.write('%s\n' % value)
