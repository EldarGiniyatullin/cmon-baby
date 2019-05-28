import numpy as np
import scipy.linalg as lg
import matplotlib.pyplot as plt


def reject_outliers(data, m=2):
    datamean = np.array([np.mean(data[:,0]), np.mean(data[:,1])])
    return np.asarray([dot for dot in data if lg.norm(dot - datamean) < m * np.std(data)])