import numpy as np
import random
import time

Samples = 0
TestCount = 0
Test2Count = 0

start_time = time.time()

print("Samples Condition[%] Test[%] Test1[%] Time[s]")

for Run in range(100000):

    N = 10 ** 8

    #print("N = ", str(N))

    for i in range(N):

        x1 = random.random()
        x2 = random.random()
        x3 = random.random()
        x4 = random.random()

        if x1 + x2 > x3 + x4:
            Samples += 1

            TestCount += x1 + x2 + x1 * x2 > x3 + x4 + x3 * x4

            Test2Count += x1 + x2 - (x3 + x4) <= x3 * x4 - x1 * x2

    #print("Number of Samples processed: ", f"{float(Samples):.2}")
    #print("Probability that the initial condition was met: ", str(Samples * 100 / (N*(Run+1))), " %")
    #print("Probability that the test was met under the intital condition: ", str(TestCount * 100 / Samples), " %")
    #print("Probability that the test2 was met under the intital condition: : ", str(Test2Count * 100 / Samples), " %")
    #print("This run took ", (time.time() - start_time), " seconds\n")

    print(f"{float(Samples):.2}", str(Samples * 100 / (N*(Run+1))), str(TestCount * 100 / Samples), str(Test2Count * 100 / Samples), f"{float(time.time() - start_time):.2}" )