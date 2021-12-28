import copy
import sys

import nose.tools as nt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.special import logsumexp
from scipy.stats import norm
from scipy.cluster import hierarchy

from viterbi import emissionLogProb, initialLogProb, transitionLogProb


def forward(df_norm, cell):
    # Positions
    bins = sorted(list(df_norm.index))

    # Set of states (copy numbers)
    Q = range(max_copy_number+1)
    
    # Initialization f[copy_number][bin] = 0
    f = [ { bin : 0 for bin in bins } for c in range(max_copy_number + 1) ]
    
    for idx, bin in enumerate(bins):
        norm_count = float(df_norm.loc[bin][cell])
        norm_count = float(df_norm.loc[bin][cell])
        for q in Q:
          eprob = emissionLogProb(q, norm_count)
          if idx == 0:
            f[q][bin] = initialLogProb(q) + eprob
          else:
            summation = []
            for copy in Q:
              summation.append(f[copy][bins[idx-1]] + transitionLogProb(copy, q))
            f[q][bin] = logsumexp(summation) + eprob

    return f

def marginal_log_prob(df_norm, f):
    # Positions
    bins = sorted(list(df_norm.index))
    
    # Set of states (copy numbers)
    Q = range(max_copy_number+1)
    last_bin = bins[-1]
    
    marginal_log_probability = None
    marginal_log_probability = logsumexp([f[c][last_bin] for c in Q])

    return marginal_log_probability