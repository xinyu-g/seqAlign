import copy
import sys

import nose.tools as nt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.special import logsumexp
from scipy.stats import norm
from scipy.cluster import hierarchy

def emissionLogProb(copy_number, norm_count):
    sigma = 0.1
    mu = copy_number / 2.

    # Compute Pr(|X_c - norm_count| <= 0.01 | c)  
    low = norm.cdf(norm_count - 0.01, mu, sigma)
    up = norm.cdf(norm_count + 0.01, mu, sigma)
    
    # prevent probability of 0
    prob = max(up - low, 0.0001)
    
    return np.log(prob)

def transitionLogProb(current_copy_number, next_copy_number):
    stay_prob = 0.99999
    if current_copy_number == next_copy_number:
        return np.log(stay_prob)
    elif 0 <= next_copy_number <= max_copy_number:
        return np.log((1 - stay_prob) / max_copy_number)
    else:
        return np.log(0)

def initialLogProb(copy_number):
    if 0 <= copy_number <= max_copy_number:
        return np.log(1./11)
    else:
        return np.log(0)


def viterbi(df_norm, cell):
    # Positions
    bins = sorted(list(df_norm.index))

    # Set of states (copy numbers)
    Q = range(max_copy_number+1)
    
    # Initialization v[copy_number][bin] = 0
    v = [ { bin : 0 for bin in bins } for c in range(max_copy_number + 1) ]
    
    for idx, bin in enumerate(bins):
        norm_count = float(df_norm.loc[bin][cell])
        for curr_copy in Q:
          eprob = emissionLogProb(curr_copy, norm_count)
          if idx == 0:
            v[curr_copy][bin] = initialLogProb(curr_copy) + eprob
          else:
            curr_max = float('-inf')
            for last_copy in Q:
              curr_max = max(curr_max, v[last_copy][bins[idx-1]] + transitionLogProb(last_copy, curr_copy) + eprob)
            v[curr_copy][bin] = curr_max

        
    return v

def max_joint_prob(df_norm, v):
    # Positions
    bins = sorted(list(df_norm.index))
    
    # Set of states (copy numbers)
    Q = range(max_copy_number+1)
    last_bin = bins[-1]
    max_joint_probability = None
    max_joint_probability = float('-inf')
    for copy in Q:
      max_joint_probability = max(max_joint_probability, v[copy][last_bin])
    

    
    return max_joint_probability

# with backtrace
def viterbi_bt(df_norm, cell):
    # Positions
    bins = sorted(list(df_norm.index))

    # Set of states (copy numbers)
    Q = range(max_copy_number+1)
    
    # Initialization v[copy_number][bin] = 0
    v = [ { bin : 0 for bin in bins } for c in range(max_copy_number + 1) ]
    bt = [ { bin : None for bin in bins } for c in range(max_copy_number + 1) ]
    
    for idx, bin in enumerate(bins):
        norm_count = float(df_norm.loc[bin][cell])
        for curr_copy in Q:
          eprob = emissionLogProb(curr_copy, norm_count)
          if idx == 0:
            v[curr_copy][bin] = initialLogProb(curr_copy) + eprob
          else:
            curr_max = float('-inf') 
            for last_copy in Q:
              p = v[last_copy][bins[idx-1]] + transitionLogProb(last_copy, curr_copy) + eprob
              if p > curr_max:
                curr_max = p
                bt[curr_copy][bin] = last_copy
            v[curr_copy][bin] = curr_max
        
    return v, bt