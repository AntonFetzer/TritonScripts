import numpy as np
from Dependencies.ReadMultipleRoot import readMultipleRoot


def EstimateError(Data, N):
    SampleLen = int(len(Data) / N)
    Sums = np.zeros(N)

    for i in range(N):
        Sums[i] = sum(Data[SampleLen * i:SampleLen * (i + 1) - 1])
        # print("From:", SampleLen*i)
        # print("To:", SampleLen*(i+1)-1)
        # print("Sum:", Samples[i])

    return np.std(Sums)*np.sqrt(N)


if __name__ == "__main__":
    Path = "/home/anton/Desktop/triton_work/6U/6U-FR4-Solder/6u-fr4-solder-2e9electrons500kev.root"
    TestData = readMultipleRoot(Path)[0]
    #TestData = TestData[0:10000]
    #print(TestData)
    print("NUmber of Data Points:", len(TestData))
    TotalDose = np.sum(TestData)
    print("TotalDose:", TotalDose)
    for x in [2, 5, 10, 20, 50, 100, 200]:
        TotalDoseStd = EstimateError(TestData, x)
        print("TotalDoseStd:", x, TotalDoseStd)

    print("sqrt(N)*Std(Data) = ", np.sqrt(len(TestData))*np.std(TestData))

    print("relative error ", 100*(np.sqrt(len(TestData))*np.std(TestData))/TotalDose, "%")

