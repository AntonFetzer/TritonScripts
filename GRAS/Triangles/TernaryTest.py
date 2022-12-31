import matplotlib.pyplot as plt
import mpltern  # noqa: F401
import numpy as np


ax = plt.subplot(projection='ternary')

Triangles = []

N = 30
Offset = (N+1)*N/2

#for x in range(N):

#t = [0, 0.5, 0.0]   # top
#l = [1, 0.5, 0.5]   # left
#r = [0, 0.0, 0.5]   # right

#ax.fill(t, l, r)

ax.set_tlabel('Top')
ax.set_llabel('Left')
ax.set_rlabel('Right')
ax.taxis.set_label_position('tick1')
ax.laxis.set_label_position('tick1')
ax.raxis.set_label_position('tick1')


h = np.sqrt(3) / 2

for x in range(N):
    for y in range(N - x):
        n = N - 1
        r = x / n
        g = y / n
        b = 1 - r - g

        if b < 0:
            print(x, y, r, g, b)
            b = 0

        ID = int(x*N+y-(x*(x-1)/2))
        #print(ID, i, ID-i, Data[ID], Colors[ID])
       # i = i + 1

        ax.fill([r+1/n, r, r], [b, b+1/n, b], [g, g, g+1/n], color=[r, g, b], linewidth=0)
        #print(x, y, r, b, g, [r+1/n, r, r], [b, b+1/n, b], [g, g, g+1/n])

for x in range(N - 1):
    for y in range(N - x - 1):
        n = N - 1
        r = x / n + 1/(3*n)
        g = y / n + 1/(3*n)
        b = 1 - r - g

        h = 1/(3*n)     # Don't ask why
        H = 2/(3*n)     # This works. I don't know why !

        if b < 0:
            print(x, y, r, g, b)
            b = 0

        ID = int(Offset+x*(N-1)+y-(x*(x-1)/2))
        #print(ID, i, ID-i, Data[ID], Colors[ID])
        #i = i + 1
        ##          top axis                Left                    Right
        #ax.fill([r+0.5/n, r+0.5/n, r], [b, b+0.5/n, b+0.5/n], [g+0.5/n, g, g+0.5/n], color=[r, g, b], linewidth=0)
        ax.fill([r-h, r+H, r+H], [b+H, b-h, b+H], [g+H, g+H, g-h], color=[r, g, b], linewidth=0)
        #print(x, y, r, b, g, [r-h, r+H, r+H], [b+H, b-h, b+H], [g+H, g+H, g-h])

plt.show()
