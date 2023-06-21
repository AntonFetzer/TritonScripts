import matplotlib.pyplot as plt
import mpltern  # noqa: F401
import numpy as np
import matplotlib as mpl
from uncertainties import ufloat
import sigfig

Fig = 0

N = 30
gamma = 0.5
Offset = int((N + 1) * N / 2)

fig = plt.figure(Fig)
ax = fig.add_subplot(projection='ternary', ternary_scale=100)

for x in range(N):
    for y in range(N - x):
        n = N - 1
        r = x / n
        g = y / n
        b = 1 - r - g

        if b < 0:
            print(x, y, r, g, b)
            b = 0

        ID = int(x * N + y - (x * (x - 1) / 2))

        r_gamma, g_gamma, b_gamma = r ** gamma, g ** gamma, b ** gamma
        max_val = max(r_gamma, g_gamma, b_gamma)
        r_normalized, g_normalized, b_normalized = r_gamma / max_val, g_gamma / max_val, b_gamma / max_val
        ax.fill([r + 1 / n, r, r], [b, b + 1 / n, b], [g, g, g + 1 / n], color=[r_normalized, g_normalized, b_normalized], linewidth=0)

for x in range(N - 1):
    for y in range(N - x - 1):
        n = N - 1
        r = x / n + 1 / (3 * n)
        g = y / n + 1 / (3 * n)
        b = 1 - r - g

        h = 1 / (3 * n)  # Don't ask why
        H = 2 / (3 * n)  # This works. I don't know why !

        if b < 0:
            print(x, y, r, g, b)
            b = 0

        ID = int(Offset + x * (N - 1) + y - (x * (x - 1) / 2))

        r_gamma, g_gamma, b_gamma = r ** gamma, g ** gamma, b ** gamma
        max_val = max(r_gamma, g_gamma, b_gamma)
        r_normalized, g_normalized, b_normalized = r_gamma / max_val, g_gamma / max_val, b_gamma / max_val
        ##          top axis                Left                    Right
        ax.fill([r - h, r + H, r + H], [b + H, b - h, b + H], [g + H, g + H, g - h], color=[r_normalized, g_normalized, b_normalized], linewidth=0)

#ax.set_title("Mass allocation map in the ternary plane", pad=20)
ax.set_tlabel("Layer 1 mass ratio [%]", color=(1, 0, 0))
ax.set_rlabel("Layer 2 mass ratio [%]", color=(0, 1, 0))
ax.set_llabel("Layer 3 mass ratio [%]", color=(0, 0, 1))

ax.taxis.set_label_position('tick1')
ax.laxis.set_label_position('tick1')
ax.raxis.set_label_position('tick1')

ax.tick_params(labelrotation='horizontal')

ax.taxis.set_tick_params(colors=(1, 0, 0), which='both')
ax.laxis.set_tick_params(colors=(0, 0, 1), which='both')
ax.raxis.set_tick_params(colors=(0, 1, 0), which='both')

for label in ax.taxis.get_ticklabels():
    label.set_color((1, 0, 0))

for label in ax.laxis.get_ticklabels():
    label.set_color((0, 0, 1))

for label in ax.raxis.get_ticklabels():
    label.set_color((0, 1, 0))

#plt.show()
plt.savefig("/home/anton/Desktop/TritonPlots/Paper/TernaryRGB.pdf", format='pdf')
