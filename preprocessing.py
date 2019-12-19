import numpy as np

def create_mask(M, start = 1, end = 99):
    p = np.sum(M, axis=1)
    mask = (p > np.percentile(p, start)) & (p < np.percentile(p, end))
    mask2D = np.outer(mask, mask)
    return mask, mask2D


def correction(M, maxiter = 1000, mask=None):
    W = M
    B = np.ones((len(M), 1))
    if mask is None:
        p = np.sum(M, axis=1)
        mask = p != 0
    converged = False
    for i in range(maxiter):
        S = np.sum(W, axis=1, keepdims=True)
        del_B = S / np.mean(S, keepdims=True)
        del_B[~mask] = 1
        del_B -= 1
        del_B *= 0.8
        del_B +=1
        mat_B = np.outer(del_B, del_B)
        W = W / mat_B
        B = B * del_B
        if np.var(S[mask]) < 1e-6:
            converged = True
            break
    corr = np.mean(B[mask])
    W = W * corr * corr
    B /= corr
    return W, B, {'converged':converged, 'iter':i}


def acd(M):
    P = np.zeros(len(M))
    for k in range(len(M)):
        P[k] = np.mean(np.diagonal(M, k))
    return P


def expected(M, P=None):
    if P is None:
        P = acd(M)
    F = np.arange(len(M))
    T = np.abs(F[None, :] - F[:, None])
    O = P[T]
    return O

def enrichment(M, O = None, mask = None):
    if O is None:
        O = expected(M)
    if mask is None:
        p = np.sum(M, axis=1)
        mask = p != 0
    mask2D = np.outer(mask, mask)
    zeros = O == 0
    W = M.copy()
    E = O.copy()

    W[zeros] = 1
    E[zeros] = 1
    enrich = W/E

    enrich = enrich[mask2D].reshape(sum(mask), sum(mask))

    return enrich
