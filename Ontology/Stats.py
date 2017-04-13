# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
Module containing functions used for statistical reasoning about the ontology
data.
"""
from Bio._py3k import range

from math import log, exp

_g = 7
_p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
     771.32342877765313, -176.61502916214059, 12.507343278686905,
     -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]


def lngamma(z):
    """
    Lanchos approximation of log((z-1)!)
    
    Reference: http://en.wikipedia.org/wiki/Lanczos_approximation
    """
    z -= 1
    x = _p[0]
    for i in range(1, _g+2):
        x += _p[i]/(z+i)
    t = z + _g + 0.5
    return 0.9189385332046727 + (z + 0.5) * log(t) - t + log(x) # log(sqrt(2*pi)) = 0.9189385332046727

def lnfactorial(n):
    """
    Returns logarithm of factorial of n.
    """
    return 0 if n < 1 else lngamma(n + 1)


def lncombination(n, k):
    """
    Returns logarithm of number of k-element sets that can be drawn from
    n-element set without replacement.
    """
    return lnfactorial(n) - lnfactorial(k) - lnfactorial(n - k)


def hypergeometric_probability(k, n, K, N):
    """
    Returns probability of k successes in n draws without replacement from
    a finite population of size N containing a maximum of K successes.
    """
    return exp(lncombination(K, k) + lncombination(N - K, n - k) - lncombination(N, n))


def hypergeometric_two_tail_test(k, n, K, N):
    """
    Returns probability of getting k or less probable number of successes in
    hypergeometric probability distribution.
    """
    um = min(n, K)
    lm = max(0, n + K - N)
    eps = 1e-15

    if um == lm:
        return 1.0

    prob = hypergeometric_probability(k, n, K, N)
    
    two_tail = 0
    for i in range(lm, um + 1):
        p = hypergeometric_probability(i, n, K, N)
        if p < prob + eps:
            two_tail += p

    two_tail = two_tail if two_tail < 1.0 else 1.0
    return two_tail

def hypergeometric_test(k, n, K, N):
    """
    Returns probability of getting k or more successes in
    hypergeometric probability distribution.
    """
    um = min(n, K)

    one_tail = hypergeometric_probability(k, n, K, N)
    
    for i in range(k, um + 1):
        one_tail += hypergeometric_probability(i, n, K, N)

    return one_tail if one_tail < 1.0 else 1.0

def bonferroni_correction(pvals):
    """
    Bonferroni correction.
    
    Reference: http://en.wikipedia.org/wiki/Bonferroni_correction
    """
    n = len(pvals)
    return [min(x * n , 1.0) for x in pvals]


def bh_fdr_correction(pvals):
    """
    Benjamin-Hochberg FDR correction.
    Reference: http://en.wikipedia.org/wiki/False_discovery_rate
    """
    n = len(pvals)
    k = 1
    cr = [1.0] * n
    mx = 0.0
    for pos, pval in sorted(enumerate(pvals), key = lambda x : x[1]):
        mx = max(pval * n / k, mx)
        if mx < 1.0:
            cr[pos] = mx
        k += 1
    return cr

def kolmogorov_smirnov_rank_test(gene_set, gene_list, gene_corr, p):
    """
    Rank test used in GSEA method. It measures dispersion of genes from
    gene_set over a gene_list. Every gene from gene_list has its weight
    specified by gene_corr. p is a parameter changing weights importance.
    
    Reference: http://www.pnas.org/content/102/43/15545.full
    """
    
    cval = 0
    Dn = 0
    Nr = 0
    
    N = len(gene_list)
    Nh = 0
    
    # Adjust correlations taking accoriding to p parameter
    if p != 1:
        adj_corr = [pow(abs(x), p) for x in gene_corr]
    else:
        adj_corr = [abs(x) for x in gene_corr]
    
    for i in range(N):
        if gene_list[i] in gene_set:
            Nr += adj_corr[i] 
            Nh += 1
    
    if N == Nh:
        miss_pen = 1.
    else:
        miss_pen = float(1) / (N - Nh)
    
    stat_plot = N * [None]
    
    for i in range(N):
        if gene_list[i] in gene_set:
            cval += adj_corr[i] / Nr
        else:
            cval -= miss_pen
        stat_plot[i] = cval
        if abs(cval) > abs(Dn):
            Dn = cval
    return (Dn, stat_plot)


def mean(arr):
    l = float(len(arr))
    return sum(arr) / l if l > .0 else .0

corrections = { "bonferroni" : bonferroni_correction,
                "bh_fdr" : bh_fdr_correction }

corrections_labels = { "bonferroni" : "Bonferroni correction",
                       "bh_fdr" : "FDR correction",
                       "fdr" : "FDR correction" }
