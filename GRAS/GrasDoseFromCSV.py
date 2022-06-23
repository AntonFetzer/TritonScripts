import matplotlib.pyplot as plt
import numpy as np
import os

Path = "/home/anton/Desktop/triton_work/GRAS-1Mat/Test/"

# Get list of all csv files in Path
CSVFiles = [f for f in os.listdir(Path + "/Results/") if f.endswith('.csv')]

print(CSVFiles)

Content = []  # 2D List to store the whole CSV file
i = 0
for file in CSVFiles:
    print(Path + Folder + "/CSV/" + file)
    Content.append([])
    with open(Path + Folder + "/CSV/" + file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            Content[i].append(row)  # Dump the whole CSV in the 2D list
    i += 1