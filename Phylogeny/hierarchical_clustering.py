import copy
import sys

import nose.tools as nt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.special import logsumexp
from scipy.stats import norm
from scipy.cluster import hierarchy


def manhattan(C, cell_1, cell_2, bins):
    value = 0

    for bin in bins:
      if C[cell_1][bin] != C[cell_2][bin]:
        value += abs(C[cell_1][bin]-C[cell_2][bin])
    return value

def cluster(distances):
    distances = copy.deepcopy(distances)
    clusters = set(distances.keys())
    n = len(clusters)
    cluster2idx = { cell : idx for idx, cell in enumerate(clusters) }
    idx2cluster = [ cell for cell in clusters ]
    Z = np.empty((0, 4), float)
    membership = [ set([cluster2idx[cell]]) for cell in clusters ]
    while len(clusters) > 1:
        # identify pair (c1, c2) with minimum distance dist
        dist, c1, c2 = None, None, None
        dist = float('inf')
        for cell_1 in clusters:
          for cell_2 in clusters:
            if cell_2 != cell_1:
              d = distances[cell_1][cell_2]
              if d < dist:
                  dist = d
                  c1 = cell_1
                  c2 = cell_2
              elif cell_1 + cell_2 < c1 + c2 and d == dist:
                  dist = d
                  c1 = cell_1
                  c2 = cell_2

        idx_c1 = cluster2idx[c1]
        idx_c2 = cluster2idx[c2]
        new_cluster = str(n + len(Z))
        new_cluster_idx = len(idx2cluster)

        membership.append( membership[idx_c1] | membership[idx_c2])
        Z = np.append(Z, 
                      np.array([[cluster2idx[c1], cluster2idx[c2], dist, len(membership[new_cluster_idx])]]), 
                      axis=0)
        
        
        distances[new_cluster] = {}
        for cell in distances.keys():
          if cell != new_cluster:
            distances[new_cluster][cell] = min(distances[cell][c1], distances[cell][c2])
            distances[cell][new_cluster] = distances[new_cluster][cell]
          else:
            distances[new_cluster][cell] = 0
            distances[cell][new_cluster] = distances[new_cluster][cell]

        clusters.add(new_cluster)
        cluster2idx[new_cluster] = len(idx2cluster)
        idx2cluster.append(new_cluster)
        
        clusters.remove(c1)
        clusters.remove(c2)
        
    return Z, idx2cluster[:n]


def complete_cluster(distances):
    distances = copy.deepcopy(distances)
    clusters = set(distances.keys())
    n = len(clusters)
    cluster2idx = { cell : idx for idx, cell in enumerate(clusters) }
    idx2cluster = [ cell for cell in clusters ]
    Z = np.empty((0, 4), float)
    membership = [ set([cluster2idx[cell]]) for cell in clusters ]
    while len(clusters) > 1:
        # identify pair (c1, c2) with minimum distance dist
        dist, c1, c2 = None, None, None

        dist = float('inf')
        for cell_1 in clusters:
          for cell_2 in clusters:
            if cell_2 != cell_1:
              d = distances[cell_1][cell_2]
              if d < dist:
                  dist = d
                  c1 = cell_1
                  c2 = cell_2
              elif cell_1 + cell_2 < c1 + c2 and d == dist:
                  dist = d
                  c1 = cell_1
                  c2 = cell_2


        idx_c1 = cluster2idx[c1]
        idx_c2 = cluster2idx[c2]
        new_cluster = str(n + len(Z))
        new_cluster_idx = len(idx2cluster)

        membership.append( membership[idx_c1] | membership[idx_c2])
        Z = np.append(Z, 
                      np.array([[cluster2idx[c1], cluster2idx[c2], dist, len(membership[new_cluster_idx])]]), 
                      axis=0)

        distances[new_cluster] = {}
        for cell in distances.keys():
          if cell != new_cluster:
            distances[new_cluster][cell] = max(distances[cell][c1], distances[cell][c2])
            distances[cell][new_cluster] = distances[new_cluster][cell]
          else:
            distances[new_cluster][cell] = 0
            distances[cell][new_cluster] = distances[new_cluster][cell]


        clusters.add(new_cluster)
        cluster2idx[new_cluster] = len(idx2cluster)
        idx2cluster.append(new_cluster)
        
        clusters.remove(c1)
        clusters.remove(c2)
      
    return Z, idx2cluster[:n]