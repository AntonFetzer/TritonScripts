import uproot
import awkward as ak
import numpy as np
import os
import multiprocessing

# PathList = ["/scratch/work/fetzera1/Gradient/2MaterialGradient/Pb-Al/root/"]
PathList = ["/home/anton/Desktop/triton_work/Gradient/3Material/Al-Pe-Cu/root/", "/home/anton/Desktop/triton_work/Gradient/3Material/Pe-Cu-Al/root/"]


def sum_up_channels(i: int, tree_, keys_: str, que: multiprocessing.Queue):
    # Edep[i] = np.sum(tree[keys[i]].array(library="np"))
    res = ak.sum(tree_[keys_[i]].array())
    print(keys_[i], res)
    que.put((i, res), timeout=1.0)


def the_job():
    for Path in PathList:

        Files = [f for f in os.listdir(Path) if f.endswith('.root')]
        # Files = ["gradient-al-test-directions-1e7proton.root"]

        for File in Files:
            process_a_file(Path, File)


def process_a_file(Path, File):
    f = uproot.open(Path + File)
    print("Read in: " + Path + File)

    tree = f["Detector Data 0"]

    keys = tree.keys()
    print(type(keys))
    print(type(keys[0]))

    # print(keys)

    Edep = np.zeros(len(keys))
    mgr = multiprocessing.Manager()
    que = mgr.Queue(1000)

    print("N cores:", multiprocessing.cpu_count())
    ncores = multiprocessing.cpu_count()-1
    nkeys = len(keys)
    processes = []
    ncjobs = min(nkeys, ncores)
    ngotten = 0
    for i in range(ncjobs):
        p = multiprocessing.Process(target=sum_up_channels, args=(i, tree, keys, que))
        p.start()
        processes.append(p)

    idx = ncjobs
    while True:
        res = que.get(timeout=50)
        ngotten += 1
        Edep[res[0]] = res[1]
        if idx < nkeys:
            p = multiprocessing.Process(target=sum_up_channels, args=(idx, tree, keys, que))
            p.start()
            processes.append(p)
            idx += 1

        if ngotten == nkeys:
            break

    for p in processes:
        p.join()

    print("Finished!")
    for i in range(len(Edep)):
        print(Edep[i])

    np.savetxt(Path + File.split(".")[0] + ".txt", Edep)

    # srun --mem=16G --time=05:00:00 python 100TilesMeVtoCSV.py


if __name__ == '__main__':
    the_job()
