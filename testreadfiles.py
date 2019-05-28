from romerolesgo import *
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import linalg as ln
from numpy import random as rd
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from outlierdetector import *

params = np.load('parameters.npz')
iterArray = array4SPS(params['u'], params['theta'], params['y'])
listOfValues = []

for i in range(144):
    npzfile = np.load('SPSLogs\\point_No_' + str(i) + '.npz')
    iterSPS = SPSmodule(iterArray, 20, 400, 25)
    iterSPS.A = npzfile['A']
    iterSPS.H = npzfile['H']
    iterSPS.g = npzfile['g']
    theta = npzfile['theta']

    listOfValues.append((theta.copy(), iterSPS.isThetaInRegionIter(theta, 1)))
    print(i)


listofDotsNonProven = [point for point, bit in listOfValues if not bit]
listOfDotsProven = [point for point, bit in listOfValues if bit]

x1, y1 = zip(*listOfDotsProven)
x2, y2 = zip(*listofDotsNonProven) # if len(listofDotsNonProven) > 0 else [], []

plt.plot(x1, y1, '.c ', params['theta'][0], params['theta'][1], 'xk ')
plt.plot(x2, y2, '.r ')

npzfile = np.load('dots_20_400_v5.npz')
print(npzfile.files)
listOdDots = npzfile['arr_0']
listOdDots1 = npzfile['arr_1']
listOdDots = reject_outliers(listOdDots)

hull = ConvexHull(listOdDots)
hull1 = ConvexHull(listOdDots1)

#plt.plot(listOdDots[:,0], listOdDots[:,1], 'o')
for simplex in hull.simplices:
    plt.plot(listOdDots[simplex, 0], listOdDots[simplex, 1], 'b-')
#plt.plot(listOdDots1[:,0], listOdDots1[:,1], 'o')
for simplex in hull1.simplices:
    plt.plot(listOdDots1[simplex, 0], listOdDots1[simplex, 1], 'g-')
plt.show()