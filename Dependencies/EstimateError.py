import numpy as np
from ..Dependencies.ReadMultipleRoot import readMultipleRoot


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
    Path = "/home/anton/Desktop/triton_work/3D/PE-W-Vault/1-5gcm2/root/pe-w-30-2e9electrons500kev.root"
    TestData = readMultipleRoot(Path)[0]
    TotalDose = np.sum(TestData)
    print("TotalDose:", TotalDose)
    for x in [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]:
        TotalDoseStd = EstimateError(TestData, x)
        print("TotalDoseStd:", x, TotalDoseStd)
