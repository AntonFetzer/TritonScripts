import numpy as np
import random
import time
import matplotlib.pyplot as plt

def f1(x, y):
    return x + y


def f2(x, y):
    return x + y + x * y


Samples = 1
TestCount = 0

start_time = time.time()

print("Samples Test[%] Diff[%] Time[s]")

for Run in range(100):

    PrevSamples = Samples
    PrevTestCount = TestCount

    N = 10 ** 8

    # print("N = ", str(N))

    for i in range(N):

        x1 = random.random()
        y1 = random.random()
        x2 = random.random()
        y2 = random.random()

        Samples += 1

        # print(x1, y1)
        # plt.scatter(x1, y1, c='b')
        # plt.scatter(x2, y2, c='r')
        if f1(x2, y2) > f1(x1, y1) and f2(x2, y2) < f2(x1, y1):
            TestCount += 1

    # print("Number of Samples processed: ", f"{float(Samples):.2}")
    # print("Probability that the initial condition was met: ", str(Samples * 100 / (N*(Run+1))), " %")
    # print("Probability that the test was met under the intital condition: ", str(TestCount * 100 / Samples), " %")
    # print("Probability that the test2 was met under the intital condition: : ", str(Test2Count * 100 / Samples), " %")
    # print("This run took ", (time.time() - start_time), " seconds\n")
    print(f"{float(Samples):.2}", str(TestCount * 100 / Samples), f"{(TestCount/Samples) - (PrevTestCount/PrevSamples):.2}", f"{float(time.time() - start_time):.2}")

#plt.show()
