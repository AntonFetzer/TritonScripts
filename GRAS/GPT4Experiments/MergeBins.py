import numpy as np

def merge_bins(hist, merge_count):
    lowerID = 0
    upperID = 1
    meanID = 2
    valueID = 3
    errorID = 4
    entriesID = 5

    if merge_count <= 1:
        return hist

    new_histogram = []
    for i in range(0, len(hist), merge_count):
        merged_bin = [0] * len(hist[0])

        for j in range(merge_count):
            if i + j < len(hist):
                merged_bin[lowerID] = hist[i][lowerID]
                merged_bin[upperID] = hist[min(i + merge_count - 1, len(hist) - 1)][upperID]
                merged_bin[meanID] += hist[i + j][meanID] * hist[i + j][entriesID]
                merged_bin[valueID] += hist[i + j][valueID]
                merged_bin[errorID] += hist[i + j][errorID] ** 2
                merged_bin[entriesID] += hist[i + j][entriesID]

        if merged_bin[entriesID] > 0:
            merged_bin[meanID] /= merged_bin[entriesID]
        merged_bin[errorID] = np.sqrt(merged_bin[errorID])

        new_histogram.append(merged_bin)

    return np.array(new_histogram)
