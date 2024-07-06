import os, pathlib, pandas, numpy, re
import scipy.sparse as sp_sparse
import scanpy as sc
from glob import glob

from pickle import FALSE, TRUE
import CytoSig
import sys

from scipy import stats
from statsmodels.stats.multitest import multipletests
from scipy.cluster import hierarchy as hc



output_path = '../out_score/cytosig/'
dataset="panE"
alpha = 1E4
nrand = 1000
alternative='two-sided'
verbose_flag = False


signature = os.path.join(sys.prefix, 'bin', 'signature.centroid')
signature = pandas.read_csv(signature, sep='\t', index_col=0)


annData = sc.read_h5ad('../data/panE.cancertype.balance.h5ad')
metadata = annData.obs
metadata = metadata.astype(object)
annData = annData.raw.to_adata()


def filter_matrix(matrix):
    
    matrix = matrix.loc[~matrix.index.isnull(), matrix.sum() > 0]
    matrix = matrix.loc[(matrix == 0).mean(axis=1) < 1]
    cnt_map = matrix.index.value_counts()
    if cnt_map.max() > 1:
        matrix = matrix.loc[cnt_map.loc[matrix.index] == 1]
    
    return matrix


def CytoSig_run_panTMEall(aData):
    

    matrix = aData.X.transpose(copy=True)
    print("matrix transpose done")
    
    matrix = pandas.DataFrame.sparse.from_spmatrix(matrix)
    print("panda dataframe done")
    
    barcodes = aData.obs_names
    genes = aData.var.index
    matrix.index = genes
    matrix.columns = barcodes
    print("matrix done")


    global metadata
    metadata=metadata.where(metadata.notnull(),'unknown')
    treatment = metadata['treatment'].apply(lambda v: v.replace('+', '_')).apply(lambda v: v.replace(' ', ''))
    info = metadata['cellType.major'] + '|' + metadata['cellType.sub'] + '|' + metadata['cancerType'] + '|'+ metadata['tissue'] + '|' + treatment + '|' + metadata['treatment_response'] + '|' + metadata['patientID'] + '|' + metadata['sampleID'].apply(lambda v: v.replace('.', '-'))
    print("colname done")

    data = filter_matrix(matrix)
    print("filter done")
    
    data.columns = info.loc[data.columns] + '|' + data.columns
    data = data.sparse.to_dense()

    
    output = os.path.join(output_path, dataset + '.signal')
    output_beta = os.path.join(output_path, dataset + '.Coef')
    output_pvalue = os.path.join(output_path, dataset + '.Pvalue')
    output_se = os.path.join(output_path, dataset + '.StdErr')
        

    background = data.mean(axis=1)
    data = data.subtract(background, axis=0)
    print("background normalize done")
    print(data.shape)
    
    beta, se, zscore, pvalue = CytoSig.ridge_significance_test(signature, data, alpha, alternative, nrand, 1, True, False, verbose_flag)
    zscore.to_csv(output, sep='\t', index_label=False)
    beta.to_csv(output_beta, sep='\t', index_label=False)
    pvalue.to_csv(output_pvalue, sep='\t', index_label=False)
    se.to_csv(output_se, sep='\t', index_label=False)
    print("CytoSig done")


CytoSig_run_panTMEall(aData=annData)






import pandas as pd
import numpy as np


ori = pd.read_csv('../out_score/cytosig/panE.signal',sep='\t')

ori.columns = annData.obs_names.tolist()
ori=ori.T

annData.obsm['cytosig_ori']=ori

sub = list(set(annData.obs['cellType.sub']))

cyto1 = pd.DataFrame(index=sub,columns=annData.obsm['cytosig_ori'].columns.tolist())

for eachc in sub:
    data = annData[annData.obs['cellType.sub']==eachc]
    for eachcy in cyto1.columns.tolist():
        cyto1.loc[eachc,eachcy]=np.mean(data.obsm['cytosig_ori'][eachcy])
        
for each in cyto1.columns.tolist():
    cyto1[each]=cyto1[each].astype('float')

cyto1.to_csv('../out_score/cytosig/celltype_cyto.csv')
