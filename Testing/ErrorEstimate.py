import numpy as np
from Dependencies.ReadMultipleRoot import readMultipleRoot
from Dependencies.ReadG4root import readG4root
import matplotlib.pyplot as plt

Path = "/home/anton/Desktop/triton_work/MULASS/MulasTrapProton2e9AlFull/root/2.root"

Data = readG4root(Path)
#Data = readMultipleRoot(Path)

Data = Data[2]

StdN = []
MaxNum = 200

for N in range(2, MaxNum):
    NumPoints = len(Data)
    SampleLen = int(NumPoints / N)

    Samples = np.zeros(N)

    for i in range(N):
        Samples[i] = sum(Data[SampleLen * i:SampleLen * (i + 1) - 1])
        # print("From:", SampleLen*i)
        # print("To:", SampleLen*(i+1)-1)
        # print("Sum:", Samples[i])

    StdN.append(np.std(Samples))


print(StdN)
x = np.arange(2, MaxNum)

plt.ylim([4, 5e3])
plt.plot(x, StdN*x, '^', label='$Std(Sums) \cdot N$')
plt.plot(x, StdN*np.sqrt(x), 'o', label='$Std(Sums) \cdot \sqrt{N}$')
plt.plot(x, StdN, 'v', label='$Std(Sums)$')
plt.yscale("log")
plt.grid(which='both')
plt.title("Convergence of error estimates")
plt.xlabel("Numer of divisions N")
plt.ylabel("Error estimate [MeV]")
plt.legend(loc='upper left')
plt.show()
