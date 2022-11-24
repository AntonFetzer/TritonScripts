import numpy as np
import matplotlib.pyplot as plt
from mpltern.ternary.datasets import get_shanon_entropies

t, l, r, v = get_shanon_entropies()

#t = [0, 0.5, 0.0]  # top
#l = [1, 0.5, 0.5]  # left
#r = [0, 0.0, 0.5]  # right
#v = [0, 1, 2]

for i in range(10):
    print(t[i], l[i], r[i], v[i])


fig = plt.figure(1)
ax = fig.add_subplot(projection='ternary')
cs = ax.tripcolor(t, l, r, v, shading='flat')
ax.set_title("TEST", pad=25)

ax.set_tlabel('Top')
ax.set_llabel('Left')
ax.set_rlabel('Right')
ax.taxis.set_label_position('tick1')
ax.laxis.set_label_position('tick1')
ax.raxis.set_label_position('tick1')

cax = ax.inset_axes([1.05, 0.1, 0.05, 0.9], transform=ax.transAxes)
colorbar = fig.colorbar(cs, cax=cax)
colorbar.set_label('Total Ionizing Dose [krad]', rotation=270, va='baseline')

plt.show()
