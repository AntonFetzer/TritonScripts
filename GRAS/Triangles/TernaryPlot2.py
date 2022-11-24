import matplotlib.pyplot as plt
import mpltern  # noqa: F401
import numpy as np


ax = plt.subplot(111, projection='ternary')

Triangles = []

N = 3

for x in range(N):

    t = [0, 0.5, 0.0]   # top
    l = [1, 0.5, 0.5]   # left
    r = [0, 0.0, 0.5]   # right

    ax.fill(t, l, r)

ax.set_tlabel('Top')
ax.set_llabel('Left')
ax.set_rlabel('Right')
ax.taxis.set_label_position('tick1')
ax.laxis.set_label_position('tick1')
ax.raxis.set_label_position('tick1')

'''

h = np.sqrt(3) / 2

for x in range(N):
    for y in range(N - x):
        n = N - 1
        r = x / n
        g = y / n
        b = 1 - r - g

        if b < 0:
            #print(x, y, r, g, b)
            b = 0

        ax.fill([r, g, b], [r + 1/N, g + 1/N, b], [r, g, b], color=([r, g, b]))
        #ax.add_patch(Polygon([(x + y / 2, y * h), (x + 1 + y / 2, y * h), (x + 0.5 + y / 2, (y + 1) * h)], color=([r, g, b])))
        #ax.text(x + y / 2 + 0.1, y * h + 0.1, str(r) + " " + str(g) + " " + str(b))
'''
plt.show()
