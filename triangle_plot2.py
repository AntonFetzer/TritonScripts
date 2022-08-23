import math
import matplotlib.pyplot as plt
import matplotlib as mpl

# big triangle max mid point
ay = 500
# big triangle center on x
ax = 0
# big triangle size edge length
Size = 1000
# big triangle height
height = math.sqrt(pow(Size,2)-pow((Size/2),2))

# sectioning as counted along one edge of the triangle
# https://stackoverflow.com/questions/65544295/creating-an-equilateral-triangle-grid-mesh-inside-a-larger-equilateral-triangl
L = 10
# triangle verticies
tri = []

dx = 0.5 * Size / L
dy = - math.sqrt(3) * dx

# generate internal triangle sectioning
for i in range(L):
    basex = ax - dx * i
    basey = ay + dy * i
    # 3 triangle points + frist one as duplicate for plotting
    tri.append([(basex, basey), (basex - dx, basey + dy), (basex + dx, basey + dy), (basex, basey)])
    for j in range(i):
        tri[-1].extend([(basex + j * 2 * dx, basey),
                        (basex + j * 2 * dx + dx, basey +   dy),
                        (basex + (j + 1) * 2 * dx, basey),
                        (basex + j * 2 * dx, basey)])
        tri[-1].extend([(basex + (j + 1) * 2 * dx, basey),
                        (basex + (j + 1) * 2 * dx - dx, basey + dy),
                        (basex + (j + 1) * 2 * dx + dx, basey + dy),
                        (basex + (j + 1) * 2 * dx, basey)])

# figure for single tri sections
fig, axis = plt.subplots(ncols=L)
# total sectioned triangle
figf, axf = plt.subplots()

for i in range(L):
    print("List len ", len(tri[i]))
    print(".")
    # get all x and y coord individually so that line plot works
    xs = []
    ys = []
    for r in range(0, len(tri[i]), 4):
        # for further sectioning visulaisation get each triangle coords
        tri_xs = []
        tri_ys = []
        for c in range(4):
            #print(c+r)
            index = c + r  
            xs.append(tri[i][index][0])
            ys.append(tri[i][index][1]) 
            tri_xs.append(tri[i][index][0])
            tri_ys.append(tri[i][index][1]) 

        # set common x, y limits
        axis[i].plot(tri_xs, tri_ys)
        axis[i].set_ylim([-(height+50)+ay, 50+ay])
        axis[i].set_xlim([-(Size/2+50)+ax, Size/2+50+ax])  


    axf.set_ylim([-(height+50)+ay, 50+ay])
    axf.set_xlim([-(Size/2+50)+ax, Size/2+50+ax])
    axf.plot(xs, ys, '-o')

plt.show()