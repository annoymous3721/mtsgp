import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from tslearn.clustering import TimeSeriesKMeans
from tslearn.datasets import CachedDatasets
from tslearn.preprocessing import TimeSeriesScalerMeanVariance, TimeSeriesResampler

seed = 0
numpy.random.seed(seed)
# Read the CSV file that generated in Rscript
df = pd.read_csv('cluster.csv')
target = pd.read_csv('target.csv')
# dtw results
X = np.array(df).reshape((9, 50, 1))
target = np.array(target).reshape((1, 50, 1))

from sdtw import SoftDTW
from sdtw.distance import SquaredEuclidean

dist_soft = np.random.rand(9)
for i in range(9):
    D = SquaredEuclidean(X[i,:,:], target[0,:,:])
    sdtw = SoftDTW(D, gamma=1.0)
    # soft-DTW discrepancy, approaches DTW as gamma -> 0
    value = sdtw.compute()
    dist_soft[i] = value

sorted_indices = np.argsort(dist_soft)
sorted_indices

# Run soft-DTW - divergence
from sdtw_div.numba_ops import sdtw_div, sdtw_div_value_and_grad

dist = np.random.rand(9)
for i in range(9):
    value = sdtw_div(X[i,:,:], target[0,:,:], gamma=1.0)
    dist[i] = value

sorted_indices = np.argsort(dist)
sorted_indices