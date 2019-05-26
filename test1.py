from romerolesgo import *
from math import sqrt
import matplotlib.pyplot as plt

def getTheRegion(h, listOfKPos, listOfKNeg, testSPS, realTheta, times, movement):
    listOfValues = []
    listOfValues10 = []
    listOfValuesTest = []
    for k in listOfKPos:
        mult = np.sqrt(1/(1+k*k))
        theta = realTheta.copy()
        thetaOut = realTheta.copy()
        thetaOut10 = realTheta.copy()
        thetaOutTest = realTheta.copy()
        isInRegion = True
        isInRegion10 = True
        isInRegionTest = True
        i = 0
        while ((isInRegion or isInRegionTest or isInRegion10) and (i < times)):
            testSPS.begin = 0
            x = theta[0] + h*mult
            y = k*(x-realTheta[0]) + realTheta[1]
            theta = np.asarray([x, y])
            if (isInRegion):
                isInRegion = True if testSPS.isThetaInRegion(theta) else False
                thetaOut = theta.copy()
            for j in range(movement):
                testSPS.isThetaInRegionHPrev1(theta)
            isInRegionTest = True if testSPS.isThetaInRegionHPrev1(theta) else False
            if (isInRegionTest):
                thetaOutTest = theta.copy()
            if (isInRegion10):
                isInRegion10 = True if testSPS.isThetaInRegion(theta) else False
                thetaOut10 = theta.copy()
            i += 1
        print("\t", i)
        listOfValues.append((thetaOut.copy(), i))
        listOfValues10.append((thetaOut10.copy(), i))
        listOfValuesTest.append((thetaOutTest.copy(), i))
        print(k)
    for k in listOfKNeg:
        mult = np.sqrt(1/(1+k*k))
        theta = realTheta.copy()
        thetaOut = realTheta.copy()
        thetaOut10 = realTheta.copy()
        thetaOutTest = realTheta.copy()
        isInRegion = True
        isInRegion10 = True
        isInRegionTest = True
        i = 0
        while ((isInRegion or isInRegionTest or isInRegion10) and (i < times)):
            testSPS.begin = 0
            x = theta[0] - h*mult
            y = k*(x-realTheta[0]) + realTheta[1]
            theta = np.asarray([x, y])
            if (isInRegion):
                isInRegion = True if testSPS.isThetaInRegion(theta) else False
                thetaOut = theta.copy()
            for j in range(movement):
                testSPS.isThetaInRegionHPrev1(theta)
            isInRegionTest = True if testSPS.isThetaInRegionHPrev1(theta) else False
            if (isInRegionTest):
                thetaOutTest = theta.copy()
            if (isInRegion10):
                isInRegion10 = True if testSPS.isThetaInRegion(theta) else False
                thetaOut10 = theta.copy()
            i += 1
        print("\t", i)
        listOfValues.append((thetaOut.copy(), i))
        listOfValues10.append((thetaOut10.copy(), i))
        listOfValuesTest.append((thetaOutTest.copy(), i))
        print(k)
    return listOfValues, listOfValues10, listOfValuesTest



listPos = np.asarray(np.arange(-30.0, -10.0, 2.0).tolist() + np.arange(-10.0, -2.0, 1.0).tolist() + np.arange(-2, 2, 0.1).tolist() + np.arange(2.0, 10.0, 1.0).tolist() + np.arange(10.0, 30.0, 2.0).tolist())
listNeg = listPos[::-1]


darealtheta = np.asarray([0.7, 0.3])
u = [np.transpose(np.asarray([0, 0]))]
y = []
for i in range(1, 45):
    u.append(np.transpose(np.asarray( [0.75 * u[i-1][0] + rd.normal(0.0, 1.0), u[i-1][0]] ) ) )
for i in range(45):
    y.append(np.dot(u[i], darealtheta) + rd.laplace(0, 0.1))

testArray1 = array4SPS(u, darealtheta, y) # rd.laplace, 0, 0.1)
testSPS1 = SPSmodule(testArray1, 20, 400, 25)
testSPS1.countRn()
print("darealtheta is", "" if testSPS1.isThetaInRegion(darealtheta) else "not", "in da region")

values, values10, valuesTest = getTheRegion(0.001, listPos, listNeg, testSPS1, darealtheta, 10000, 9)

listOfDots = [point for point, bit in values if bit]
listOfDots10 = [point for point, bit in values10 if bit]
listOfDotsTest = [point for point, bit in valuesTest if bit]

x1, y1 = zip(*listOfDots)
x2, y2 = zip(*listOfDots10)
x3, y3 = zip(*listOfDotsTest)

plt.plot(x1, y1, '.b ', x2, y2, '.g ', x3, y3, '.r ', darealtheta[0], darealtheta[1], 'xk ')

plt.show()
F = open("test_results_matrix.txt", "w")
for t in listOfDots:
  F.write(' '.join(str(s) for s in t) + '\n')
F.write("\n")
for t in listOfDots10:
  F.write(' '.join(str(s) for s in t) + '\n')
F.write("\n")
for t in listOfDotsTest:
  F.write(' '.join(str(s) for s in t) + '\n')
F.write("\n")
F.close()
plt.savefig('C:\\Space\\Progs\\PishiDavai\\foo.png')

