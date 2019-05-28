from romerolesgo import *
from math import sqrt
import matplotlib.pyplot as plt

def getTheRegion(_h, listOfKPos, listOfKNeg, testSPS, realTheta, times, makeLog=False, logfile=None):
    listOfValues = []
    iter = 0
    for k in listOfKPos:
        mult = np.sqrt(1/(1+k*k))
        theta = realTheta.copy()
        isInRegion = True
        i = 0
        #TODO: automatic system here
        while (isInRegion and (i < times)):
            testSPS.begin = 0
            x = theta[0] + _h*mult
            y = k*(x-realTheta[0]) + realTheta[1]
            theta = np.asarray([x, y])
            isInRegion = True if testSPS.isThetaInRegion(theta) else False
            i += 1
        if makeLog:
            np.savez(logfile + str(iter) + '.npz', theta=theta, A=testSPS.A, R=testSPS.Rn, H=testSPS.H, g=testSPS.g)
        iter += 1
        print("\t", i)
        listOfValues.append((theta.copy(), i))
        print(k)
    for k in listOfKNeg:
        mult = np.sqrt(1/(1+k*k))
        theta = realTheta.copy()
        isInRegion = True
        i = 0
        while (isInRegion and (i < times)):
            testSPS.begin = 0
            x = theta[0] - _h*mult
            y = k*(x-realTheta[0]) + realTheta[1]
            theta = np.asarray([x, y])
            isInRegion = True if testSPS.isThetaInRegion(theta) else False
            i += 1
        if makeLog:
            np.savez(logfile + str(iter) + '.npz', theta=theta, A=testSPS.A, R=testSPS.Rn, H=testSPS.H, g=testSPS.g)
            iter += 1
        print("\t", i)
        listOfValues.append((theta.copy(), i))
        print(k)
    return listOfValues



listPos = np.asarray(np.arange(-30.0, -6.0, 2.0).tolist() + np.arange(-6.0, -2.0, 1).tolist() +
                     np.arange(-2, 2, 0.05).tolist() + np.arange(2.0, 6.0, 1).tolist() +
                     np.arange(6.0, 30.0, 2.0).tolist())
listNeg = listPos[::-1]


darealtheta = np.asarray([0.7, 0.3])
u = [np.transpose(np.asarray([0, 0]))]
y = []
for i in range(1, 45):
    u.append(np.asarray( [0.75 * u[i-1][0] + rd.normal(0.0, 1.0), u[i-1][0]] )  )
for i in range(45):
    y.append(np.dot(u[i], darealtheta) + rd.laplace(0, 0.1))
np.savez('parameters.npz', u=u, y=y, theta=darealtheta)

testArray1 = array4SPS(u, darealtheta, y) # rd.laplace, 0, 0.1)
testSPS1 = SPSmodule(testArray1, 20, 400, 25)
testSPS1.countRn()
print("darealtheta is", "" if testSPS1.isThetaInRegion(darealtheta) else "not", "in da region")

h = 0.0005
stepCount = 3000
values = getTheRegion(h, listPos, listNeg, testSPS1, darealtheta, stepCount, False)#True, 'SPSLogs\\point_No_')

listOfDots = [point for point, i in values if i < stepCount]
listofDotsOverCount = [point for point, i in values if i >= stepCount]

x1, y1 = zip(*listOfDots)
x2, y2 = zip(*listofDotsOverCount) if listofDotsOverCount else [], []

plt.plot(x1, y1, '.b ', darealtheta[0], darealtheta[1], 'xk ')
plt.plot(x2, y2, '.r ')

testSPS1.begin += 20
values20 = getTheRegion(h, listPos, listNeg, testSPS1, darealtheta, stepCount)

listOfDots2 = [point for point, i in values20 if i < stepCount]
listofDotsOverCount2 = [point for point, i in values20 if i >= stepCount]

x1_2, y1_2 = zip(*listOfDots2)
x2_2, y2_2 = zip(*listofDotsOverCount2) if listofDotsOverCount2 else [], []

plt.plot(x1_2, y1_2, '.g ')
plt.plot(x2_2, y2_2, '.k ')

np.savez('dots_20_400_v8.npz', listOfDots, listOfDots2)

plt.show()

