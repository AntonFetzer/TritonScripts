from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
from GRAS.Dependencies.TotalKRadGras import totalkRadGras
from matplotlib import cm

Path = "/home/anton/Desktop/triton_work/3MatTriangles/PE-Al-Pb/Res/"

#Electrons = totalkRadGras(Path, "Elec")
Protons = totalkRadGras(Path, "Prot")
#Total = Electrons + Protons
Total = Protons

ColorData = Total[0]
Min = np.min(ColorData)
ColorData = ColorData - np.min(ColorData)
Max = np.max(ColorData)
ColorData = ColorData/Max

Colors = cm.turbo(ColorData)
print(len(ColorData))

N = 30
Offset = int((N+1)*N/2)

print("The average Dose is", np.mean(Total[0]))
print("The average Dose of the main triangle is", np.mean(Total[0][:Offset]))
print("The average Dose of the upside down triangle is", np.mean(Total[0][Offset:]))
print("The difference in Dose is", 100*(np.mean(Total[0][Offset:]) - np.mean(Total[0][:Offset])) / np.mean(Total[0]), "%")

h = np.sqrt(3) / 2

fig, ax = plt.subplots(figsize=(20, 18))

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

        ax.add_patch(Polygon([(x + y / 2, y * h), (x + 1 + y / 2, y * h), (x + 0.5 + y / 2, (y + 1) * h)], color=Colors[ID], linewidth=0))
        #ax.text(x + y / 2 + 0.1, y * h + 0.1, str(r) + " " + str(g) + " " + str(b))


for x in range(N - 1):
    for y in range(N - x - 1):
        n = N - 1
        r = x / n + 1/(3*n)
        g = y / n + 1/(3*n)
        b = 1 - r - g

        if b < 0:
            #print(x, y, r, g, b)
            b = 0

        ID = int(Offset+x*(N-1)+y-(x*(x-1)/2))
        #print(ID, i, ID-i, Data[ID], Colors[ID])
        #i = i + 1

        ax.add_patch(Polygon([(x + y / 2 + 0.5, (y + 1) * h), (x + 1 + y / 2, y * h), (x + 1.5 + y / 2, (y + 1) * h)], color=Colors[ID], linewidth=0))
        #ax.text(x + y / 2 + 0.5, (y + 1) * h - 0.2, str(r) + " " + str(g) + " " + str(b))


plt.ylim(0, N * h)
plt.xlim(0, N)
plt.axis('off')
plt.gca().set_aspect('equal')
#plt.colorbar()
plt.title("The average Dose is " + str(np.mean(Total[0])) + " krad")

#plt.savefig(Path + "../TIDmap.eps", format='eps', bbox_inches="tight")
plt.show()
