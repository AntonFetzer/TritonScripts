import numpy as np

Data = np.zeros(32*32)

A = 1
B = 2
C = 3

for x in range(32):
    for y in range(100):
        Data[x+100*y] = A*(31-x)*(31-y) + B*x**(31-y) +

print(Data)