import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import linalg as ln
from numpy import random as rd
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from outlierdetector import *

npzfile = np.load('dots_20_400_v7_manyKs.npz')
print(npzfile.files)
listOdDots = npzfile['arr_0']
listOdDots1 = npzfile['arr_1']
listOdDots = reject_outliers(listOdDots)
listOdDots1 = reject_outliers(listOdDots1)

hull = ConvexHull(listOdDots)
hull1 = ConvexHull(listOdDots1)

plt.plot(listOdDots[:,0], listOdDots[:,1], 'o')
for simplex in hull.simplices:
    plt.plot(listOdDots[simplex, 0], listOdDots[simplex, 1], 'b-')
plt.plot(listOdDots1[:,0], listOdDots1[:,1], 'o')
for simplex in hull1.simplices:
    plt.plot(listOdDots1[simplex, 0], listOdDots1[simplex, 1], 'g-')
plt.show()