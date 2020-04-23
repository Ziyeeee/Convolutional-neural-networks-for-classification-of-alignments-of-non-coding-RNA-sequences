import numpy as np
import sys

import time

def assembledata(ddir):

    checkdata = np.load(ddir+"/portion/ncRNApair_data0.npy")
    dsize = len(checkdata) + 1
    genelen = len(checkdata[0])
    width = len(checkdata[0][0])

    data = np.empty((0,genelen,width), dtype=np.float32)
    label = np.empty((0), dtype=np.int32)

    print("Makeing data...")
    start = time.time()

    for i in range(dsize):
        indata = ddir + "/portion/ncRNApair_data"+str(i)+".npy"
        inlabel = ddir + "/portion/ncRNApair_label"+str(i)+".npy"

        tmpdata = np.load(indata)
        tmplabel = np.load(inlabel)
        
        print(tmpdata.shape[0])

        data = np.concatenate([data, tmpdata], axis=0)
        label = np.concatenate([label, tmplabel], axis=0)

    print("Complete makeing data.")
    elapsed_time = time.time() - start
    print ("Loding time:{0}".format(elapsed_time) + "[sec]")
    print("Data size :",len(data),"*",len(data[0]),"*",len(data[0][0]))

    np.save(ddir+"/ncRNApair_data.npy", data)
    np.save(ddir+"/ncRNApair_labe.npy", label)


if __name__ == '__main__':

    datadir = sys.argv[1]
    if datadir[-1] == "/":
        datadir = datadir[:-1]

    assembledata(datadir)
