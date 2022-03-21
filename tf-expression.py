import pandas as pd
import numpy as np
from functools import reduce
from sklearn.linear_model import LassoCV
from sklearn.metrics import r2_score
from collections import defaultdict
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import fdrcorrection
from random import sample
from sklearn.linear_model import LinearRegression
from argparse import ArgumentParser

# read the names of input and output files
parser = ArgumentParser(description='Get parameters')
parser.add_argument('-g', dest='gFile', type=str,
                   help='Normalized gene expression values')
parser.add_argument('-p', dest='pFile', type=str,
                   help='Result from Predixcan')
parser.add_argument('-o', dest='oFile', type=str,
                   help='Output file')

if len(sys.argv) != 3:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

gene_expression_path = ReadParams(args.gFile)

# remove covariants
cov = pd.read_table(f'{tissue}.v8.covariates.txt', index_col=0)
covariates = cov.iloc[list(range(3))+list(range(5, 20))+[-2, -1],:]


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
l = []
for i, chunk in enumerate(expression):
    new = remove_covariate_component(chunk, covariates)
    l.append(new)
removed = pd.concat(l).reset_index()

# Cov-removed gene expression and the result of predixcan
prediction_path = ReadParams(args.pFile)
prediction = pd.read_table(prediction_path)
cols =[x.split('.')[0] for x in list(prediction.columns)]
prediction.columns = cols
predixcan = prediction.set_index('IID').T
removed = removed.set_index('gene_id')
common_col = list(set((removed.columns).intersection(predixcan.columns)))
removed1 = removed[common_col]
predixcan1 = predixcan[common_col]
common_gene = list(set((removed1.index).intersection(predixcan1.index)))
removed2 = removed1[removed1.index.isin(common_gene)]
predixcan2 = predixcan1[predixcan1.index.isin(common_gene)]
predixcan2 = predixcan2.loc[removed2.index, :]

# Calculate the R2 scores of predixcan genes
# Select genes with R2 scores equal to or large than 0.1
r2_scores = {}
for i, row in removed2.iterrows():
    y = row.values.reshape(-1, 1)
    y_p = predixcan2.loc[i, :].values.reshape(-1, 1)
    lm = LinearRegression()
    lm.fit(y_p, y)
    y_new = lm.predict(y_p)
    r2 = 1-((y_new-y)**2).sum()/((y-y.mean())**2).sum()
    r2_scores[i] = r2
r2_score_csv = pd.DataFrame(r2_scores.items(), columns=['gene', 'r2_score'])
r2_scores = dict(zip(r2_score_csv['gene'], r2_score_csv['r2_score']))
genes = list(dict(filter(lambda x: x[1] >=0.1, r2_scores.items())).keys())

normalized = removed2
predixcan = predixcan2
diff = normalized - predixcan

# Robustness check
'''
all_samples = diff.columns
sample_length = len(all_samples)
random_length = int(0.9 * sample_length)
random_col = sample(list(all_samples), random_length)
diff = diff[random_col]
'''

expression = pd.read_csv(gene_expression_path, sep='\t')
expression.gene_id = expression.gene_id.apply(lambda x: x.split('.')[0])
expression = expression.iloc[:, 3:]
expression = expression.set_index('gene_id')
expression = expression[diff.columns]

# Mapping between gene codes and gene name
ENSG_to_name_table = pd.read_csv('ENSG_to_name.txt')
ENSG_to_gene_name = dict(zip(ENSG_to_name_table['Gene stable ID'], ENSG_to_name_table['HGNC symbol']))
gene_name_to_ENSG = dict(zip(ENSG_to_name_table['HGNC symbol'], ENSG_to_name_table['Gene stable ID']))

# Mapping between genes and their TFs
tf2gene = pd.read_table('TF2Gene.txt')
tf2gene.TF = tf2gene.TF.apply(lambda x: gene_name_to_ENSG.get(x, None))
tf2gene.Gene = tf2gene.Gene.apply(lambda x: gene_name_to_ENSG.get(x, None))
ENSG_to_tfs = dict(tf2gene.groupby('Gene').apply(lambda x: list(x.TF)))


# Filter out the TFs without expression values
def remove_tf(tf_list):
    new_tfs = []
    for tf in tf_list:
        if tf in expression.index:
            new_tfs.append(tf)
    return new_tfs


ENSG_to_tfs = {k: remove_tf(v) for k, v in ENSG_to_tfs.items()}
ENSG_to_tfs = {k: v for k, v in ENSG_to_tfs.items() if len(v)>0}


def build_model(gene, expression, diff, ENSG_to_tfs):
    if gene not in ENSG_to_tfs:
        #print(f'{gene} does not have TFs')
        return None, None
    tfs = ENSG_to_tfs[gene]
    X = expression.loc[tfs, :].T
    y = diff.loc[gene, :]
    lassocv = LassoCV(cv=10, random_state=1)
    lassocv.fit(X,y)
    y_pred = lassocv.predict(X)
    r2 = r2_score(y, y_pred)
    return lassocv, r2


models = {}
r2s = {}
for gene in genes:
    models[gene], r2s[gene] = build_model(gene, expression, diff, ENSG_to_tfs)
    
    
# Shuffle the target
def shuffle_y(gene, expression, diff, ENSG_to_tfs):
    shuffle_100_r2 = []
    if gene not in ENSG_to_tfs:
        return None
    tfs = ENSG_to_tfs[gene]
    X = expression.loc[tfs, :].T
    y = diff.loc[gene, :]
    for i in range(100):
        y_shuffle = np.random.permutation(y)
        while any(y_shuffle == y):
            y_shuffle = np.random.permutation(y)
        lassocv = LassoCV(cv=10, random_state=1)
        lassocv.fit(X,y_shuffle)
        y_pred = lassocv.predict(X)
        r2 = r2_score(y_shuffle, y_pred)
        shuffle_100_r2.append(r2)
    return shuffle_100_r2


p_values = defaultdict(float)
for gene in genes:
    gene_r2_list = shuffle_y(gene, expression, diff, ENSG_to_tfs)
    if isinstance(gene_r2_list, list):
        # ztest
        _, ztest_pvalue = ztest(x1=np.array([r2s[gene]]), x2=np.array(gene_r2_list), alternative='larger')
        p_values[gene] = ztest_pvalue

# Calculate the adjusted p values.        
_, p_values_corrected = fdrcorrection(np.asarray(list(p_values.values())))
p_values_new = dict(zip(p_values.keys(), p_values_corrected))
genes_after_shuffle_y = {k:v for k, v in sorted(p_values_new.items(), key=lambda x: x[0]) if v < 0.05}
genes_after_shuffle_y_list = list(genes_after_shuffle_y.keys())


# Shuffle the TFs of the gene
def generate_shuffled_tf(all_tfs_set, ENSG_to_tfs):
    shuffled_ENSG_to_tfs = defaultdict(list)
    for gene, tfs in ENSG_to_tfs.items():
        shuffled_ENSG_to_tfs[gene] = sample(all_tfs_set.difference(tfs), len(tfs))
    return shuffled_ENSG_to_tfs


def shuffle_tf(gene, expression, diff, ENSG_to_tfs):
    all_tfs_set = reduce(lambda a, b: set(a).union(b), list(ENSG_to_tfs.values()))
    shuffle_100_r2 = []
    for i in range(100):
        shuffled_ENSG_to_tfs = generate_shuffled_tf(all_tfs_set, ENSG_to_tfs)
        _, r2 = build_model(gene, expression, diff, shuffled_ENSG_to_tfs)
        shuffle_100_r2.append(r2)
    return shuffle_100_r2


p_values_tf = defaultdict(float)
for gene in genes_after_shuffle_y_list:
    gene_r2_list = shuffle_tf(gene, expression, diff, ENSG_to_tfs)
    if isinstance(gene_r2_list, list):
        # ztest
        _, ztest_pvalue = ztest(x1=np.array([r2s[gene]]), x2=np.array(gene_r2_list), alternative='larger')
        p_values_tf[gene] = ztest_pvalue
tmp, p_values_corrected_tf = fdrcorrection(np.asarray(list(p_values_tf.values())))
p_values_new_tf = dict(zip(p_values_tf.keys(), p_values_corrected_tf))
genes_after_shuffle_tf = {k:v for k, v in sorted(p_values_new_tf.items(), key=lambda x: x[0]) if v < 0.05}

# Result after background models
significant_genes = [ENSG_to_gene_name[x] for x in list(genes_after_shuffle_tf)]
print('Gene | R2 score of predixcan | R2 score of our model')
for ensg in genes_after_shuffle_tf.keys():
    print(f'{ENSG_to_gene_name[ensg]}: {r2_scores[ensg]} | {r2s[ensg]}')  
tf_for_sig = {}
for g in list(genes_after_shuffle_tf.keys()):
    tfs = [ENSG_to_gene_name[i] for i in ENSG_to_tfs[g]]
    tf_for_sig[ENSG_to_gene_name[g]] = sorted(tfs)
for k,v in tf_for_sig.items():
    print(k ,v )
tf_weights = {}
for i in list(genes_after_shuffle_tf.keys()):
    tfs = [ENSG_to_gene_name[j] for j in ENSG_to_tfs[i]]
    weights = models[i].coef_
    weights_dict = dict(zip(tfs, weights))
    tf_weights[ENSG_to_gene_name[i]] = weights_dict 
for k, v in tf_weights.items():
    with open(ReadParams(args.oFile), 'a+') as f:
        f.write(f'=====Gene: {k}=====')
        f.write('TF and Weights: ')
        f.write(v)