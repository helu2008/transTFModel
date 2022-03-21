import pandas as pd
import numpy as np
from sklearn.linear_model import LassoCV
from sklearn.metrics import r2_score
from matplotlib import pyplot as plt
from collections import defaultdict
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import fdrcorrection
from functools import reduce
from random import sample
from pandas_plink import read_plink1_bin
from sklearn.linear_model import LinearRegression

# read the names of input and output files
parser = ArgumentParser(description='Get parameters')
parser.add_argument('-g', dest='gFile', type=str,
                   help='Normalized gene expression values')
parser.add_argument('-p', dest='pFile', type=str,
                   help='Result from Predixcan')
parser.add_argument('-bed', dest='bedFile', type=str,
                   help='genotype bed file')
parser.add_argument('-bim', dest='bimFile', type=str,
                   help='genotype bim file')
parser.add_argument('-fam', dest='famFile', type=str,
                   help='genotype fam file')
parser.add_argument('-snpSift', dest='snpSiftFile', type=str,
                   help='SNPs after filtered by SnpSift')
parser.add_argument('-s', dest='sFile', type=str,
                   help='SNPs after filtered by Sift Score')
parser.add_argument('-o', dest='oFile', type=str,
                   help='Output file')

if len(sys.argv) != 8:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

# Gene expression value
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

# Calculate the R2 scores and select genes with R2 scores equal to or large than 0.1
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
genes_with_2_75 = [x for x in genes if len(ENSG_to_tf_variant.get(x,set()))<=75 and len(ENSG_to_tf_variant.get(x,set()))>=2]
genes = genes_with_2_75

normalized = removed2
predixcan = predixcan2
diff = normalized - predixcan
diff = diff.reset_index()
expression = pd.read_csv(gene_expression_path, sep='\t')
expression.gene_id = expression.gene_id.apply(lambda x: x.split('.')[0])
expression = expression.iloc[:, 3:]
expression = expression.set_index('gene_id')
expression = expression[diff.columns]

# Mapping between gene codes and gene name
ENSG_to_name_file = pd.read_csv('ENSG_to_name.txt')
ENSG_to_name = dict(zip(ENSG_to_name_file['Gene stable ID'], ENSG_to_name_file['HGNC symbol']))
name_to_ENSG = dict(zip(ENSG_to_name_file['HGNC symbol'], ENSG_to_name_file['Gene stable ID']))

# Mapping from gene to its TFs
tf2gene = pd.read_table('TF2Gene.txt')
gene2tf = dict(tf2gene.groupby('Gene').apply(lambda x: list(x.TF)))
tf2gene.TF = tf2gene.TF.apply(lambda x: gene_name_to_ENSG.get(x, None))
tf2gene.Gene = tf2gene.Gene.apply(lambda x: gene_name_to_ENSG.get(x, None))
ENSG_to_tfs = dict(tf2gene.groupby('Gene').apply(lambda x: list(x.TF)))

# Extract genotype data
nssnp_genotype = read_plink1_bin(ReadParams(args.bedFile), ReadParams(args.bimFile), ReadParams(args.famFile), verbose=False)

# Mapping
genotype = pd.read_csv(ReadParams(args.bimFile), sep='\t', header=None)
chrom_pos_to_name = {}
for i, row in genotype.iterrows():
    chrom = row[0]
    pos = row[3]
    name = row[1]
    ref = row[4]
    alt = row[5]
    if len(ref) > 1 or len(alt) > 1:
        continue
    chrom_pos_to_name[(chrom, pos)] = name 
pos_to_variant = dict(zip(genotype[1], genotype.index))

# Create range file
with open(ReadParams(args.snpSiftFile), 'r') as f:
    rsid_chrom_pos = defaultdict(tuple)
    line = f.readline()
    while line.startswith('#'):
        line = f.readline()
    for line in f:
        info = line.strip().split('\t')
        chrom = info[0]
        pos = info[1]
        rsid = info[2]
        rsid_chrom_pos[rsid] = (chrom, pos)
with open('range_rsid.txt', 'w') as f:
    for k, v in rsid_chrom_pos.items():
        chrom = v[0]
        pos = v[1]
        rsid = k
        f.write(f'{chrom} {pos} {pos} {rsid}\n')
range_file = pd.read_csv('range_rsid.txt', sep=' ', header=None)
rsid_to_chrom_pos = {}
for i, row in range_file.iterrows():
    chrom = row[0]
    pos = row[1]
    rsid = row[3]
    if chrom == 'X':
        chrom = 23
    elif chrom == '22':
        chrom = 22
    elif chrom in {'Y', 'MT'}:
        continue
    chrom = int(chrom)
    if (chrom, pos) in chrom_pos_to_name:
        snp_name = chrom_pos_to_name[(chrom, pos)]
        rsid_to_chrom_pos[rsid] = snp_name
rsid_to_variant = {k: pos_to_variant[v] for k, v in rsid_to_chrom_pos.items()}
variant_to_rsid = {v: k for k, v in rsid_to_variant.items()}

# filter deleterious SNPs using sift score
with open(ReadParams(args.sFile), 'r') as f:
    gene_nssnp_list_del = defaultdict(set)
    line = f.readline()
    while line.startswith('#'):
        line = f.readline()
    for i, line in enumerate(f):
        nssnp_id = line.split('\t')[2]
        info = line.split('\t')[-1].split(';')
        for j in info:
            if j.startswith('GENEINFO'):
                gene = j.split('=')[1].split(':')[0]
                gene_nssnp_list_del[gene].add(nssnp_id)
                
gene_name_to_rsids_before_pph = gene_nssnp_list_del                
# input list of tfs
# return set of nssnps of the tfs
def merge_nssnps(l):
    s = set()
    for t in l:
        if t in gene_name_to_rsids_before_pph:
            tmp = gene_name_to_rsids_before_pph.get(t)
            s.update(tmp)
    return s
ENSG_tf_rsids = dict(filter(lambda x: x[0] and x[1], {gene_name_to_ENSG.get(k) : merge_nssnps(v) for k,v in gene2tf.items()}.items()))
def rsid_set_to_var_set(s):
    result = set()
    for i in s:
        v = rsid_to_variant.get(i)
        if v is not None:
            result.add(v)
    return result
ENSG_to_tf_variant = dict(filter(lambda x: x[1],{k: rsid_set_to_var_set(v) for k, v in ENSG_tf_rsids.items()}.items()))

# Build model
nssnp_model = {}
r2score_dict = {}
for gene in genes:
    if gene in ENSG_to_tf_variant:
        X = nssnp_genotype.sel(sample=list(diff.columns)[1:], variant=["variant"+str(x) for x in ENSG_to_tf_variant[gene]]).values
        X = np.where(np.isnan(X), 0, X)
        y = diff[diff.gene_id==gene].values[0, 1:].astype('d')
        lassocv = LassoCV(cv=10, random_state=0)
        lassocv.fit(X,y)
        y_pred = lassocv.predict(X)
        r2 = r2_score(y, y_pred)
        r2score_dict[ENSG_to_name[gene]] = r2
        nssnp_model[ENSG_to_name[gene]] = lassocv
 
  
def shuffle_y_build_model_10(gene):
    shuffle_list = []
    if gene in ENSG_to_tf_variant and len(ENSG_to_tf_variant[gene])>0:
        X = nssnp_genotype.sel(sample=list(diff.columns)[1:], variant=["variant"+str(x) for x in ENSG_to_tf_variant[gene]]).values
        X = np.where(np.isnan(X), 0, X)
        y = diff[diff.gene_id==gene].values[0, 1:].astype('d')
        for i in range(10):
            y_shuffle = np.random.permutation(y)
            while any(y_shuffle == y):
                y_shuffle = np.random.permutation(y)
            lassocv = LassoCV(cv=10, random_state=0)
            lassocv.fit(X,y_shuffle)
            y_pred = lassocv.predict(X)
            r2 = r2_score(y_shuffle, y_pred)
            shuffle_list.append(r2)
        return shuffle_list

# permute y
p_value = defaultdict(float)
for gene in genes:
    if gene in ENSG_to_tf_variant:
        shuffle_list = shuffle_y_build_model_10(gene)
        _, ztest_pvalue = ztest(x1=np.array([r2score_dict[ENSG_to_name[gene]]]), x2=np.array(shuffle_list), alternative='larger')
        p_value[gene] = ztest_pvalue

_, p_value_corrected = fdrcorrection(np.asarray(list(p_value.values())))
p_value_new = dict(zip(p_value.keys(), p_value_corrected))
genes_5 = {k:v for k, v in sorted(p_value_new.items(), key=lambda x: x[0]) if v < 0.05}

var_set = reduce(lambda a, b: a.union(b), ENSG_to_tf_variant.values())
def generate_shuffled_nssnp():
    shuffled_var = {}
    for gene, v in ENSG_to_tf_variant.items():
        if len(ENSG_to_tf_variant[gene])>0:
            shuffled_var[gene] = sample(var_set.difference(v), len(v))
    return shuffled_var

def shuffle_nssnp_build_model_10(gene):
    shuffle_list = []
    if gene in ENSG_to_tf_variant:
        for j in range(10):
            shuffled_var = generate_shuffled_nssnp()
            X = nssnp_genotype.sel(sample=list(diff.columns)[1:], variant=["variant"+str(x) for x in shuffled_var[gene]]).values
            X = np.where(np.isnan(X), 0, X)
            y = diff[diff.gene_id==gene].values[0, 1:].astype('d')
            lassocv = LassoCV(cv=10, random_state=0)
            lassocv.fit(X,y)
            y_pred = lassocv.predict(X)
            r2 = r2_score(y, y_pred)
            shuffle_list.append(r2)
        return shuffle_list
    
# permute nssnp
p_value2 = defaultdict(float)
for gene in genes:
    if gene in genes_5:
        if gene in ENSG_to_tf_variant and len(ENSG_to_tf_variant[gene])>0:
            shuffle_list = shuffle_nssnp_build_model_10(gene)
            _, ztest_pvalue = ztest(x1=np.array([r2score_dict[ENSG_to_name[gene]]]), x2=np.array(shuffle_list), alternative='larger')
            p_value2[gene] = ztest_pvalue

_, p_value_corrected2 = fdrcorrection(np.asarray(list(p_value2.values())))
p_value_new2 = dict(zip(p_value2.keys(), p_value_corrected2))
genes_5_nssnp = {k:v for k, v in sorted(p_value_new2.items(), key=lambda x: x[0]) if v < 0.05}

sig_genes_10_rounds = list(genes_5_nssnp.keys())
genes = sig_genes_10_rounds

def shuffle_y_build_model(gene):
    shuffle_list = []
    if gene in ENSG_to_tf_variant and len(ENSG_to_tf_variant[gene])>0:
        X = nssnp_genotype.sel(sample=list(diff.columns)[1:], variant=["variant"+str(x) for x in ENSG_to_tf_variant[gene]]).values
        X = np.where(np.isnan(X), 0, X)
        y = diff[diff.gene_id==gene].values[0, 1:].astype('d')
        for i in range(100):
            y_shuffle = np.random.permutation(y)
            while any(y_shuffle == y):
                y_shuffle = np.random.permutation(y)
            lassocv = LassoCV(cv=10, random_state=0)
            lassocv.fit(X,y_shuffle)
            y_pred = lassocv.predict(X)
            r2 = r2_score(y_shuffle, y_pred)
            shuffle_list.append(r2)
        return shuffle_list
    
# permute y
p_value = defaultdict(float)
for gene in genes:
    if gene in ENSG_to_tf_variant:
        shuffle_list = shuffle_y_build_model(gene)
        _, ztest_pvalue = ztest(x1=np.array([r2score_dict[ENSG_to_name[gene]]]), x2=np.array(shuffle_list), alternative='larger')
        p_value[gene] = ztest_pvalue
        
_, p_value_corrected = fdrcorrection(np.asarray(list(p_value.values())))
p_value_new = dict(zip(p_value.keys(), p_value_corrected))
genes_5 = {k:v for k, v in sorted(p_value_new.items(), key=lambda x: x[0]) if v < 0.05}

def shuffle_nssnp_build_model(gene):
    shuffle_list = []
    if gene in ENSG_to_tf_variant:
        for j in range(100):
            shuffled_var = generate_shuffled_nssnp()
            X = nssnp_genotype.sel(sample=list(diff.columns)[1:], variant=["variant"+str(x) for x in shuffled_var[gene]]).values
            X = np.where(np.isnan(X), 0, X)
            y = diff[diff.gene_id==gene].values[0, 1:].astype('d')
            lassocv = LassoCV(cv=10, random_state=0)
            lassocv.fit(X,y)
            y_pred = lassocv.predict(X)
            r2 = r2_score(y, y_pred)
            shuffle_list.append(r2)
        return shuffle_list
    
    
# permute nssnp
p_value2 = defaultdict(float)
for gene in genes:
    if gene in genes_5:
        if gene in ENSG_to_tf_variant and len(ENSG_to_tf_variant[gene])>0:
            shuffle_list = shuffle_nssnp_build_model(gene)
            _, ztest_pvalue = ztest(x1=np.array([r2score_dict[ENSG_to_name[gene]]]), x2=np.array(shuffle_list), alternative='larger')
            p_value2[gene] = ztest_pvalue
    
tmp, p_value_corrected2 = fdrcorrection(np.asarray(list(p_value2.values())))
p_value_new2 = dict(zip(p_value2.keys(), p_value_corrected2))
genes_5_nssnp = {k:v for k, v in sorted(p_value_new2.items(), key=lambda x: x[0]) if v < 0.05}

significant_gene_list = []
for k in genes_5_nssnp.keys():
    significant_gene_list.append(ENSG_to_name[k])

rsid_to_gene_name = dict()
for k,v in gene_name_to_rsids_before_pph.items():
    for i in v:
        rsid_to_gene_name[i]=k
variant_to_rsid_mapping = {v: k for k, v in rsid_to_variant.items()}

genes = list(genes_5_nssnp.keys())
snp_weight = {}
print('Gene | R2 score of predixcan | R2 score of our model')
for ensg in genes_5_nssnp.keys():
    print(f'{ENSG_to_gene_name[ensg]}: {r2_scores[ensg]} | {r2score_dict[ENSG_to_gene_name[ensg]]}')  
for gene in [ENSG_to_name[i] for i in genes]:
    gene2 = name_to_ENSG[gene]
    weight_dict = defaultdict(lambda:0)
    model = nssnp_model[gene]
    variants = list(ENSG_to_tf_variant[gene2])
    for i in range(len(variants)):
        tf = rsid_to_gene_name[variant_to_rsid_mapping[variants[i]]]
        weight_dict[tf] = model.coef_[i] if abs(weight_dict.get(tf, 0)) < abs(model.coef_[i]) else weight_dict.get(tf, 0)
    weight_dict = sorted(weight_dict.items(), key=lambda x: abs(x[1]), reverse=True)
    snp_weight[gene] = dict(weight_dict)
significant_gene_tf_nssnp = {}
for gene in snp_weight.keys():
    model = nssnp_model[gene]
    gene = name_to_ENSG[gene]
    variants = list(ENSG_to_tf_variant[gene])
    rsid_list = [variant_to_rsid_mapping[i] for i in variants]
    tf_rsid_tuple = [(rsid_to_gene_name[i], i) for i in rsid_list]
    rsid_weights_dict = dict(zip(tf_rsid_tuple, model.coef_))
    significant_gene_tf_nssnp[ENSG_to_name[gene]] = rsid_weights_dict
for k, v in significant_gene_tf_nssnp.items():
    with open(ReadParams(args.oFile), 'a+') as f:
        f.write(f'=====Gene: {k}=====')
        f.write('SNPs and Weights: ')
        f.write(v)