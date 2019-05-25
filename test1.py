from romerolesgo import *
from math import sqrt
import matplotlib.pyplot as plt

darealtheta = np.asarray([0.7, 0.3])
u = [np.transpose(np.asarray([0, 0]))]
y = []
for i in range(1, 35):
    u.append(np.transpose(np.asarray( [0.75 * u[i-1][0] + rd.normal(0.0, 1.0), u[i-1][0]] ) ) )
for i in range(35):
    y.append(np.dot(u[i], darealtheta) + rd.laplace(0, 0.1))

testArray1 = array4SPS(u, darealtheta, rd.laplace, 0, 0.1)
testSPS1 = SPSmodule(testArray1, 20, 400, 25)
testSPS1.countRn()
test1realTheta95 = testSPS1.isThetaInRegion(darealtheta)


'''trues = 0
trys = 3

testSPS1.countRn()
for i in range(trys):
    testSPS1.generateA()
    test1_95_2 = testSPS1.isThetaInRegion(np.asarray([1.4, 1.4]))
    if test1_95_2:
        trues += 1
print(trues, "of", trys, "are True")
'''

listPos = np.asarray(np.arange(-20, -2, 1).tolist() + np.arange(-2, 2, 0.1).tolist() + np.arange(2, 20, 1).tolist())
listNeg = np.asarray(np.arange(20, 2, -1).tolist() + np.arange(2, -2, -0.1).tolist() + np.arange(-2, -20, -1).tolist())

def getTheRegion(h, listOfKPos, listOfKNeg, testSPS, realTheta):
    listOfValues = []
    for k in listOfKPos:
        mult = np.sqrt(1/(1+k*k))
        theta = realTheta
        isInRegion = True
        i = 0
        while (isInRegion and (i < 1000)):
            x = theta[0] + h*mult
            y = k*(x-darealtheta[0]) + darealtheta[1]
            theta = np.asarray([x, y])
            test1_95_2 = testSPS.isThetaInRegion(theta)
            isInRegion = True if test1_95_2 else False
            i += 1
            print("\t", i)
        listOfValues.append((theta.copy(), i))
        print(k)
    for k in listOfKNeg:
        mult = np.sqrt(1/(1+k*k))
        theta = realTheta
        isInRegion = True
        i = 0
        while (isInRegion and (i < 1000)):
            x = theta[0] - h*mult
            y = k*(x-darealtheta[0]) + darealtheta[1]
            theta = np.asarray([x, y])
            test1_95_2 = testSPS.isThetaInRegion(theta)
            isInRegion = True if test1_95_2 else False
            i += 1
        print("\t", i)
        listOfValues.append((theta.copy(), i))
        print(k)
    return listOfValues

listOfValues = getTheRegion(0.01, listPos, listNeg, testSPS1, darealtheta)

listOfDots1 = [point for point, bit in listOfValues if bit]
listOfDots2 = [point for point, bit in listOfValues if not bit]

x1, y1 = zip(*listOfDots1)
if not (len(listOfDots2) is 0):
    x2, y2 = zip(*listOfDots2)
    plt.plot(x1, y1, '.b ', x2, y2, '.r ', darealtheta[0], darealtheta[1], 'xk ')
else:
    plt.plot(x1, y1, '.b ', darealtheta[0], darealtheta[1], 'xk ')

plt.show()
plt.savefig('C:\\Space\\Progs\\PishiDavai\\foo.png')
