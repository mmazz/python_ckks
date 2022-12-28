import numpy as np
import random

""" Diferentes distributions used in the ckkks implementation."""
arith64bits = True

def sampleRing(ql, N):
    """ Samples a vector in the ring space."""
    ql2 = np.rint(ql/2)
    coef = np.random.randint(-ql2+1,ql2+1, N, dtype=np.int64)
    return np.rint(coef)

def hammingWeight(coef: np.array):
    """ Calculate the hamming weight of a vector."""
    h = 0
    for i in range(len(coef)):
        if(coef[i]!=0):
            h+=1
    return h

def sample_Sparce_ternary(h, N):
    """ Set of signed binary vector whose Hamming weight is exactly h.
    This implementation is by brute force...
    """
    h_test = 0
    while(h_test!=h):
        coef = np.random.randint(-1, 1, N, dtype=np.int64)
        h_test = hammingWeight(coef)
    return coef

def sample_Uniform_ternary(N):
    """ Set of signed binary uniformly distributed."""
    coef = np.random.randint(-1, 2, N, dtype=np.int64)
    return coef

def sampleSecret(h, N, sparce_ternary=False):
    """ Samples a signed binary vector.
    The original paper used a sparce ternary distribution, but for security
    reason it change to a uniform ternary distribution.
    """
    if sparce_ternary:
        coef = sample_Sparce_ternary(h, N)
    else:
        coef = sample_Uniform_ternary(N)
    return coef


def sampleDG(sigma, N):
    """ Samples to a discrete gaussian of variance sigma^2."""
    coef = np.int64(np.random.normal(loc=0.0, scale=sigma, size=N))
    return np.rint(coef)


def sampleZO(N, rho=0.5):
    """ Samples a distribution that has a probability 1-rho to be zero
    and a probability rho/2 for each of  1 and -1.
    Where 0<= rho <= 1.
    In this version its hardcoded to rho=0.5 like the paper.
    """
    coef = np.zeros(N)
    for i in range(len(coef)):
        a = random.randint(0, 1)
        if a==1:
            a = random.randint(0, 1)
            if a==0:
                a = -1
            coef[i] = a
    return np.int64(coef)

