import numpy as np
import time

for i in range(10):
    start_time = time.time()
    N = 10 ** i

    print("N = ", str(N))

    x1 = np.random.rand(N)
    x2 = np.random.rand(N)
    x3 = np.random.rand(N)
    x4 = np.random.rand(N)

    Condition = x1 + x2 > x3 + x4

    Test = (x1 + x2 + x1 * x2 > x3 + x4 + x3 * x4) & Condition

    Test2 = (x1 + x2 - (x3 + x4) <= x3 * x4 - x1 * x2) & Condition

    ProbCondition = np.sum(Condition) / N
    ProbTest = np.sum(Test) / N / ProbCondition
    ProbTest2 = np.sum(Test2) / N / ProbCondition

    print("Probability that the initial condition was met: ", str(ProbCondition * 100), " %")
    print("Probability that the test condition was met: ", str(ProbTest * 100), " %")
    print("Probability that the test2 condition was met: ", str(ProbTest2 * 100), " %")
    print("This run took ", (time.time() - start_time), " seconds")

# plt.plot(Condition, '.')
# plt.plot(Test, '.')
# plt.show()
