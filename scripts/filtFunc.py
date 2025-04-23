#import
from scipy.stats import binom, binomtest
import graph_tool.all as gt
import leidenalg as la
import logging as lg
import igraph as ig
import pandas as pd
import numpy as np
logger = lg.getLogger()
logger.setLevel("DEBUG")

#scaler
def scaler(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def jaccard_similarity(g, v1, v2):
    n1 = set(g.neighbors(v1))
    n2 = set(g.neighbors(v2))
    intersection = len(n1.intersection(n2))
    union = len(n1.union(n2))
    if union == 0:
        return 0
    else:
        return intersection / union

#noise_corrector
def noise_corrector(g, sdev_score = False):
    nc = g.copy()
    if nc.is_directed():
        print("Implement directed version.")
        return
    nc.es["p_value"] = np.nan
    nc.es["sdev"] = np.nan
    nc.es["score"] = np.nan
    n = np.sum(nc.es["weight"])
    for e in nc.es:
        i = e.source
        j = e.target
        w = e["weight"]
        ni = nc.vs[i].strength(weights = "weight", mode = "out")
        nj = nc.vs[j].strength(weights = "weight", mode = "in")
        mean_prior_probability = ((ni * nj) / n) * (1 / n)
        kappa = n / (ni * nj)
        nc.es[e.index]["p_value"] = 1 - binom.cdf(w, n, mean_prior_probability)
        if sdev_score:
            score = ((kappa * w) - 1) / ((kappa * w) + 1)
            var_prior_probability = (1 / (n ** 2)) * (ni * nj * (n - ni) * (n - nj)) / ((n ** 2) * ((n - 1)))
            alpha_prior = (((mean_prior_probability ** 2) / var_prior_probability) * (1 - mean_prior_probability)) - mean_prior_probability
            beta_prior = (mean_prior_probability / var_prior_probability) * (1 - (mean_prior_probability ** 2)) - (1 - mean_prior_probability)
            alpha_post = alpha_prior + w
            beta_post = n - w + beta_prior
            expected_pij = alpha_post / (alpha_post + beta_post)
            variance_nij = expected_pij * (1 - expected_pij) * n
            d = (1.0 / (ni * nj)) - (n * ((ni + nj) / ((ni * nj) ** 2)))
            variance_cij = variance_nij * (((2 * (kappa + (w * d))) / (((kappa * w) + 1) ** 2)) ** 2)
            sdev_cij = variance_cij ** 0.5
            nc.es[e.index]["sdev"] = sdev_cij
            nc.es[e.index]["score"] = score
    return nc

#marginal_likelihood
def marginal_likelihood(g, max_neg_log = 0):
    ml = g.copy()
    if ml.is_directed():
        print("Implement directed version.")
        return
    ks = ml.strength(weights = "weight")
    td = np.sum(ks)
    for e in ml.es:
        i0, i1 = e.source, e.target
        try:
            w = e["weight"]
            ku = ks[i0]
            kv = ks[i1]
            q = td / 2.0
            p = ku * kv * 1.0 / q / q / 2.0
            p = binomtest(k = int(w), n = int(q), p = p, alternative = "greater").pvalue
            e["p_value"] = max(max_neg_log, max_neg_log if p <= 0 else p)
        except ValueError as error:
            logger.warning("warning: ValueError {}".format(str(error)))
            logger.debug("ValueError weight: {} ks[i0]: {} ks[i1]: {} td: {} p: {}".format(e["weight"], ks[i0], ks[i1], td, p))
            e["p_value"] = None
        except Exception as error:
            logger.warning("warning: Exception {}".format(str(error)))
            e["p_value"] = None
    max_sig = np.min([s for s in ml.es["p_value"] if s is not None])
    for e in ml.es:
        if e["p_value"] is None:
            e["p_value"] = max_sig
    return ml
