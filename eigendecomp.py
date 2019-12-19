import numpy as np
from scipy.sparse.linalg import eigsh
import preprocessing
import pandas as pd


def get_eig(M, n=10, stats = False):
    S, U = eigsh(M, n)
    order = np.argsort(-np.abs(S))
    S = S[order]
    U = U.T[order]
    if stats:
        print(np.round(S/np.sum(np.abs(S)), 5)[:3])
        print('First 3 eigenvectors represent',
          np.round(np.abs(S[:3]).sum()/np.sum(np.abs(S)), 5),
          'of the ' +str(n)+ ' eigenvectors')
    return S, U.T

def get_outer(M, n= 10, S=None, U = None):
    if (S is None) or (U is None):
        S, U = get_eig(M,n)
    return np.array( [S[i]*np.outer(U[:, i],U[:, i]) for i in range(len(U.T))])



def eigendecomp(M, operation='clip', n =3, genetrack = None, mask = None):
    enrich = preprocessing.enrichment(M, mask = mask)
    enrich = apply_op(enrich, operation)
    if enrich is None:
        return None
    S, U = get_eig(enrich, n)
    if not (genetrack is None):
        for i in range(n):
            U[:,i] = orient_vector(U[:,i], genetrack)
    return S, U

def cooler_eigendecomp(cool, bins, regions, operations = ['clip'], trackname = 'gene_count', n = 3):

    df_list  =list()

    for ch in regions:
        mat = cool.matrix(balance = False).fetch(ch)
        mask, mask2D = preprocessing.create_mask(mat, start = 1, end = 99)
        mat[~mask2D] = 0
        W, B, res = preprocessing.correction(mat, maxiter = 10000)
        loc_eig = bins[bins.chrom == ch].copy()
        genes = loc_eig[trackname].values[mask]
        Z = [eigendecomp(W, operation=op, n=n, genetrack = genes, mask = mask) for op in operations]
        i = 0
        for op in operations:

            for k in range(n):
                track = np.zeros(len(loc_eig))
                track[:] = np.nan
                track[mask] = Z[i][1][:,k]
                loc_eig['E' + str(k+1)+ '_' + op] = track
            i+=1
        df_list.append(loc_eig)

    return pd.concat(df_list)


def orient_vector(vector, genetrack):
    p = correlation1D( vector, genetrack)
    if p < 0:
        return -vector
    return vector

def correlation1D(vector1, vector2):
    g = vector1.reshape(-1)
    v = vector2.reshape(-1)
    dots = np.dot(v - np.mean(v), v - np.mean(v)) * np.dot(g - np.mean(g), g - np.mean(g))

    rel = np.dot(g - np.mean(g),v - np.mean(v))
    return rel/(dots**0.5)

def correlation2D(M1, M2):
    return correlation1D(M1.reshape(-1), M2.reshape(-1))

def apply_op(M, operation = 'clip'):
    if operation == 'clip':
        return np.clip(M, 0, np.percentile(M, 99)) -1
    elif operation == 'log':
        return np.log(M + 1e-1)
    elif operation == 'cliplog':
        return np.clip(np.log(M + 1e-1), np.percentile(np.log(M + 1e-1), 1), np.percentile(np.log(M + 1e-1), 99))
    elif operation == 'tanh':
        return np.tanh(M-1)
    elif operation == 'tanhlog':
        return np.tanh(np.log(M + 1e-1))
    else:
        print('Ilegal operation')
        return None
