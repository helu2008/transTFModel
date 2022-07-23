import sys
import pandas as pd
import numpy as np
from functools import reduce
from sklearn.linear_model import LassoCV
from sklearn.linear_model import LassoLarsCV
from sklearn.metrics import r2_score
from collections import defaultdict
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import fdrcorrection
from random import sample
from sklearn.linear_model import LinearRegression
from argparse import ArgumentParser
import timeit
import warnings
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning

simplefilter("ignore", category=ConvergenceWarning)

#warnings.filterwarnings('ignore')

# general params
CUTOFF_PVALUE = 0.05
NUM_SHUFFLES = 100
NUM_ROBUST_ITER = 10
FRACTION_ROBUST_SAMPLES = 0.9

# read the names of input and output files
parser = ArgumentParser(description='Get parameters')
parser.add_argument('-g', dest='gFile', type=str,
                   help='Normalized gene expression values')
parser.add_argument('-p', dest='pFile', type=str,
                   help='Result from Predixcan')
parser.add_argument('-c', dest='covFile', type=str,
                   help='Tissue covariance file')
parser.add_argument('-o', dest='oFile', type=str,
                   help='Output file')

#print (sys.argv)
if len(sys.argv) < 9:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

gene_expression_path = args.gFile

# remove covariants

start = timeit.default_timer()
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

    
full_expression = pd.read_table(gene_expression_path, chunksize=100)
chunk_residuals = []
for i, chunk in enumerate(full_expression):
    new = remove_covariate_component(chunk, covariates)
    chunk_residuals.append(new)
removed = pd.concat(chunk_residuals).reset_index()

# Cov-removed gene expression and the result of predixcan
prediction_path = args.pFile
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

end = timeit.default_timer()
print('Finished removing covariates. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

# Calculate the R2 scores of predixcan genes
# Select genes with R2 scores equal to or large than 0
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
genes = list(dict(filter(lambda x: x[1] >=0, r2_scores.items())).keys())

normalized = removed2
predixcan = predixcan2
diff = normalized - predixcan

end = timeit.default_timer()
print('Finished calculating R2. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

full_expression = pd.read_csv(gene_expression_path, sep='\t')
full_expression.gene_id = full_expression.gene_id.apply(lambda x: x.split('.')[0])
full_expression = full_expression.iloc[:, 3:]
full_expression = full_expression.set_index('gene_id')
expression = full_expression[diff.columns]

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

end = timeit.default_timer()
print('Finished mapping TF expression. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

def build_model(gene, expression, diff, ENSG_to_tfs):
    if gene not in ENSG_to_tfs:
        #print(f'{gene} does not have TFs')
        return None, None
    tfs = ENSG_to_tfs[gene]
    X = expression.loc[tfs, :].T
    y = diff.loc[gene, :]
    #lassocv = LassoLarsCV(cv=10, n_jobs=1, normalize=True, eps=1e-10)  # ,random_state=1
    lassocv = LassoCV(cv=10, n_jobs=1, tol=1e-2)  # ,random_state=1, normalize=True,
    lassocv.fit(X,y)
    y_pred = lassocv.predict(X)
    r2 = r2_score(y, y_pred)
    #lassocv2 = LassoCV(cv=10, n_jobs=1, normalize=True, eps=1e-10)  # ,random_state=1
    #lassocv2.fit(X,y)
    #y_pred2 = lassocv2.predict(X)
    #r2_2 = r2_score(y, y_pred2)
    return lassocv, r2

def build_model_new (gene, expression, diff, tfs):
    X = expression.loc[tfs, :].T
    y = diff.loc[gene, :]
    #lassocv = LassoLarsCV(cv=10, n_jobs=1, normalize=True, eps=1e-10)  # ,random_state=1
    lassocv = LassoCV(cv=10, n_jobs=1, tol=1e-2)  # ,random_state=1, normalize=True,
    lassocv.fit(X,y)
    y_pred = lassocv.predict(X)
    r2 = r2_score(y, y_pred)
    #lassocv2 = LassoCV(cv=10, n_jobs=1, normalize=True, eps=1e-10)  # ,random_state=1
    #lassocv2.fit(X,y)
    #y_pred2 = lassocv2.predict(X)
    #r2_2 = r2_score(y, y_pred2)
    return lassocv, r2

models = {}
r2s = {}
counter = 0
#f = open("compare_lasso.txt", 'w')
#f.write("Regular\tLARS")

for gene in genes:
    counter = counter + 1
#    if counter % 10 == 0:
#        end = timeit.default_timer()
#        print('Finished gene model for ' + gene + '. ' + str(end - start) + ' seconds')
#        start = timeit.default_timer()
    models[gene] = None
    r2s[gene] = None
    if gene in ENSG_to_tfs:
        models[gene], r2s[gene] = build_model_new (gene, expression, diff, ENSG_to_tfs[gene])
    #f.write(str(r2_2) + "\t" + str(r2s[gene]) + "\n")

# remove genes with model results of 0 or None
filteredGenes = dict(filter(lambda x: x[1] is not None, r2s.items()))
genes = list(dict(filter(lambda x: x[1] > 0, filteredGenes.items())))
end = timeit.default_timer()
print('Finished modeling genes. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

# Shuffle the target
def shuffle_y(gene, expression, diff, ENSG_to_tfs):
    shuffles_r2 = []
    if gene not in ENSG_to_tfs:
        return None
    tfs = ENSG_to_tfs[gene]
    X = expression.loc[tfs, :].T
    y = diff.loc[gene, :]
    for i in range(NUM_SHUFFLES):
        y_shuffle = np.random.permutation(y)
#        while any(y_shuffle == y):
#            y_shuffle = np.random.permutation(y)
        #lassocv = LassoLarsCV(cv=10, n_jobs=1, normalize=True, eps=1e-10)
        lassocv = LassoCV(cv=10, random_state=1, tol=1e-2)
        lassocv.fit(X,y_shuffle)
        y_pred = lassocv.predict(X)
        r2 = r2_score(y_shuffle, y_pred)
        shuffles_r2.append(r2)
    return shuffles_r2


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

genes_after_shuffle_y = {k:v for k, v in sorted(p_values_new.items(), key=lambda x: x[0]) if v < CUTOFF_PVALUE}
genes_after_shuffle_y_list = list(genes_after_shuffle_y.keys())

end = timeit.default_timer()
print('Finished shuffling genes. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

# Shuffle the TFs of the gene
def generate_shuffled_tf(all_tfs_set, ENSG_to_tfs):
    shuffled_ENSG_to_tfs = defaultdict(list)
    for gene, tfs in ENSG_to_tfs.items():
        shuffled_ENSG_to_tfs[gene] = sample(list(all_tfs_set.difference(tfs)), len(tfs))
    return shuffled_ENSG_to_tfs


def shuffle_tf(gene, expression, diff, ENSG_to_tfs):
    if gene not in ENSG_to_tfs:
        return None
    all_tfs_set = reduce(lambda a, b: set(a).union(b), list(ENSG_to_tfs.values()))
    tfs = ENSG_to_tfs[gene]
    shuffle_tf_r2 = []
    for i in range(NUM_SHUFFLES):

        shuffled_tfs = sample(list(all_tfs_set.difference(tfs)), len(tfs))
        # shuffled_ENSG_to_tfs = generate_shuffled_tf(all_tfs_set, ENSG_to_tfs)
        r2 = 0.0
        _, r2 = build_model_new(gene, expression, diff, shuffled_tfs)
        shuffle_tf_r2.append(r2)
    return shuffle_tf_r2


p_values_tf = defaultdict(float)
for gene in genes_after_shuffle_y_list:
    gene_r2_list = shuffle_tf(gene, expression, diff, ENSG_to_tfs)
    if isinstance(gene_r2_list, list):
        # ztest
        _, ztest_pvalue = ztest(x1=np.array([r2s[gene]]), x2=np.array(gene_r2_list), alternative='larger')
        p_values_tf[gene] = ztest_pvalue
tmp, p_values_corrected_tf = fdrcorrection(np.asarray(list(p_values_tf.values())))
p_values_new_tf = dict(zip(p_values_tf.keys(), p_values_corrected_tf))
genes_after_shuffle_tf = {k:v for k, v in sorted(p_values_new_tf.items(), key=lambda x: x[0]) if v < CUTOFF_PVALUE}

end = timeit.default_timer()
print('Finished shuffling TFs. Took ' + str(int(end-start)) + ' seconds')
start = timeit.default_timer()

# Result after background models
significant_genes = [ENSG_to_gene_name[x] for x in list(genes_after_shuffle_tf)]
with open(args.oFile, 'w') as f:
    f.write('Gene\tR2 score of predixcan\tR2 score of our model\n')
    for ensg in genes_after_shuffle_tf.keys():
        f.write(f'{ENSG_to_gene_name[ensg]}\t{r2_scores[ensg]}\t{r2s[ensg]}\n')
    #tf_for_sig = {}
    #for g in list(genes_after_shuffle_tf.keys()):
    #    tfs = [ENSG_to_gene_name[i] for i in ENSG_to_tfs[g]]
    #    tf_for_sig[ENSG_to_gene_name[g]] = sorted(tfs)
    #for k, v in tf_for_sig.items():
    #    f.write('%s\t%s\n' % (k, v))
tf_weights = {}
for i in list(genes_after_shuffle_tf.keys()):
    tfs = [ENSG_to_gene_name[j] for j in ENSG_to_tfs[i]]
    weights = models[i].coef_
    weights_dict = dict(zip(tfs, weights))
    tf_weights[ENSG_to_gene_name[i]] = weights_dict
with open(args.oFile, 'a+') as f:
    f.write('----------------------------\nTF and Weights: \n')
    f.write('Gene\tTFs\n')
    for k, v in tf_weights.items():
        f.write('%s\t%s\n' % (k, v))
        #f.write(str(v) + '\n')

# Robustness check
genes = list(genes_after_shuffle_tf.keys())
random90hitgenes = []

# expression = pd.read_csv(gene_expression_path, sep='\t')
# expression.gene_id = expression.gene_id.apply(lambda x: x.split('.')[0])
# expression = expression.iloc[:, 3:]
# expression = expression.set_index('gene_id')
# expression = expression[diff.columns]
for i in range(NUM_ROBUST_ITER):
    print ("Robustness check #" + i)
    diff = normalized - predixcan
    all_samples = diff.columns
    sample_length = len(all_samples)
    random_length = int(FRACTION_ROBUST_SAMPLES * sample_length)
    random_col = sample(list(all_samples), random_length)
    diff = diff[random_col]
    expression = full_expression[diff.columns]
    models = {}
    r2s = {}
    for gene in genes:
        # print(gene)
        if gene in ENSG_to_tfs:
            models[gene], r2s[gene] = build_model_new(gene, expression, diff, ENSG_to_tfs[gene])
        # print(r2s[gene])
        
    p_values = defaultdict(float)
    for gene in genes:
        gene_r2_list = shuffle_y(gene, expression, diff, ENSG_to_tfs)
        if isinstance(gene_r2_list, list):
            # ztest
            _, ztest_pvalue = ztest(x1=np.array([r2s[gene]]), x2=np.array(gene_r2_list), alternative='larger')
            p_values[gene] = ztest_pvalue
            #print(f'{ENSG_to_gene_name[gene]}: {ztest_pvalue}')
            
    _, p_values_corrected = fdrcorrection(np.asarray(list(p_values.values())))
    p_values_new = dict(zip(p_values.keys(), p_values_corrected))
    genes_after_shuffle_y = {k:v for k, v in sorted(p_values_new.items(), key=lambda x: x[0]) if v < CUTOFF_PVALUE}
    genes_after_shuffle_y_list = list(genes_after_shuffle_y.keys())
    
    p_values_tf = defaultdict(float)
    for gene in genes_after_shuffle_y_list:
        gene_r2_list = shuffle_tf(gene, expression, diff, ENSG_to_tfs)
        if isinstance(gene_r2_list, list):
            # ztest
            _, ztest_pvalue = ztest(x1=np.array([r2s[gene]]), x2=np.array(gene_r2_list), alternative='larger')
            p_values_tf[gene] = ztest_pvalue
            #print(f'{ENSG_to_gene_name[gene]}: {ztest_pvalue}')
            
    _, p_values_corrected_tf = fdrcorrection(np.asarray(list(p_values_tf.values())))
    p_values_new_tf = dict(zip(p_values_tf.keys(), p_values_corrected_tf))
    genes_after_shuffle_tf = {k:v for k, v in sorted(p_values_new_tf.items(), key=lambda x: x[0]) if v < CUTOFF_PVALUE}
    significant_genes = [ENSG_to_gene_name[x] for x in list(genes_after_shuffle_tf)]
    random90hitgenes.append(significant_genes)

end = timeit.default_timer()
print('Finished robustness checks. Took ' + str(int(end-start)) + ' seconds')

with open(args.oFile, 'a+') as f:
    f.write("--------------\nRobustness check\n")
    for value in random90hitgenes:
        f.write('%s\n' % value)
    # f.write(str(random90hitgenes))