import uproot
import awkward as ak
import numpy as np
import os
import multiprocessing

# PathList = ["/scratch/work/fetzera1/Gradient/2MaterialGradient/Pb-Al/root/"]
PathList = ["/scratch/work/fetzera1/Gradient/2Material1-5gcm2/Al-Al/root/"]


def sum_up_channels(i: int, tree_, keys_: str, que: multiprocessing.Queue):
    # Edep[i] = np.sum(tree[keys[i]].array(library="np"))
    res = ak.sum(tree_[keys_[i]].array())
    # print(keys_[i], res)
    que.put((i, res), timeout=10)


def sum_up_channels2(jobque: multiprocessing.Queue, resque: multiprocessing.Queue, tree, keys):
    print("\tprocess thread {} started".format(multiprocessing.Process.pid), flush=True)
    while True:
        try:
            job = jobque.get(timeout=60)
        except:
            raise TimeoutError("Sum up timeout")
        i = job
        if i == -1:
            return
        r = ak.sum(tree[keys[i]].array())
        resque.put((i, r), timeout=10)


def process_a_file(Path, File):
    mgr = multiprocessing.Manager()

    f = uproot.open(Path + File, array_cache="20000 MB")
    print("Read in: " + Path + File)
    tree = f["Detector Data 0"]
    keys = tree.keys()

    print("Number of cores =", multiprocessing.cpu_count())
    nprocs = multiprocessing.cpu_count() - 1
    nkeys = len(keys)
    processes = dict()

    Edep = np.zeros(len(keys))

    resque = mgr.Queue(nkeys + 2)
    jobque = mgr.Queue(nkeys + 2)

    print("Starting {} processes...".format(nprocs))
    for i in range(nprocs):
        p = multiprocessing.Process(target=sum_up_channels2, args=(jobque, resque, tree, keys))
        #p = threading.Thread(target=sum_up_channels2, args=(jobque, resque, tree, keys))
        p.start()
        processes[i] = p

    print("Filling jobque...")
    for i in range(nkeys):
        jobque.put(i)

    print("Gathering results...")
    ngot = 0
    while ngot < nkeys:
        (i, res) = resque.get(timeout=60)
        Edep[i] = res
        print("Got result idx:{} = {}".format(i, res))
        ngot += 1

    print("Killing threads...")
    for i in range(nprocs):
        jobque.put(-1)

    for idx, p in processes.items():
        p.join()

    f.close()

    print("Finished!")

    np.savetxt(Path + File.split(".")[0] + ".txt", Edep)

    # srun --mem=16G --time=05:00:00 python 100TilesMeVtoCSV.py


def the_job():
    for Path in PathList:

        Files = [f for f in os.listdir(Path) if f.endswith('.root')]
        # Files = ["3layer2e9electron.root"]

        for File in Files:
            process_a_file(Path, File)


if __name__ == '__main__':
    the_job()
