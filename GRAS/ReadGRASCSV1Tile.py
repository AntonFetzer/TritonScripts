import csv


def readGrasCsv1Tile(file):
    data = {}  # Initialise empty dict !

    ReadFlag = 0
    # print("Reading in File: " + file)
    with open(file, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            # print(line)
            if ReadFlag == 0:
                if "'Non zero entries'" in line:
                    ReadFlag = 1
            elif ReadFlag == 1:
                if "'End of Block'" in line:
                    ReadFlag = 2
                else:
                    data["Dose"] = float(line[0])
                    data["Error"] = float(line[1])
                    data["Entries"] = float(line[2])
                    data["NonZero"] = float(line[3])
    return data


if __name__ == "__main__":
    File = "/home/anton/Desktop/triton_work/GRAS-1Mat/ParallelTest/Results/Electrons_4325_268.csv"

    print(readGrasCsv1Tile(File))
