from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

N = 5

Offset = (N+1)*N/2

h = np.sqrt(3) / 2

fig, ax = plt.subplots()

i = 0

for x in range(N):
    for y in range(N - x):
        n = N - 1
        r = x / n
        g = y / n
        b = 1 - x / n - y / n

        if b < 0:
            #print(x, y, r, g, b)
            b = 0

        ID = int(x*N+y-(x*(x-1)/2))
        #print(ID, i, ID-i, Data[ID], Colors[ID])
        #i = i + 1

        ax.add_patch(
            Polygon([(x + y / 2, y * h), (x + 1 + y / 2, y * h), (x + 0.5 + y / 2, (y + 1) * h)], color=[r, g, b]))
        #ax.text(x + y / 2 + 0.1, y * h + 0.1, str(r) + " " + str(g) + " " + str(b))
        ax.text(x + y / 2 + 0.3, y * h + 0.2, str([r, g, b]))


for x in range(N - 1):
    for y in range(N - 1 - x):
        n = N - 1
        r = x / n + 0.5 / n
        g = y / n + 0.5 / n
        b = 1 - x / n - y / n - 1 / n

        if b < 0:
            #print(x, y, r, g, b)
            b = 0

        ID = int(Offset+x*(N-1)+y-(x*(x-1)/2))
        #print(ID, i, ID-i, Data[ID], Colors[ID])
        #i = i + 1

        ax.add_patch(Polygon([(x + y / 2 + 0.5, (y + 1) * h), (x + 1 + y / 2, y * h), (x + 1.5 + y / 2, (y + 1) * h)], color=[r, g, b]))
        #ax.text(x + y / 2 + 0.5, (y + 1) * h - 0.2, str(r) + " " + str(g) + " " + str(b))
        ax.text(x + y / 2 + 0.7, (y + 1) * h - 0.3, str([r, g, b]))


plt.ylim(0, N * h)
plt.xlim(0, N)
plt.axis('off')
plt.gca().set_aspect('equal')
#plt.colorbar()

plt.show()
