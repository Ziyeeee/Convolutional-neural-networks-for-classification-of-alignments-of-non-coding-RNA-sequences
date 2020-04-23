# coding=utf-8

import argparse
import math

import numpy as np
import chainer
from chainer.datasets import tuple_dataset

import time
import gc

randomseed = 777

# listをn個のsublistに等分割する 
def split_list (list,n):
    list_size = len(list)
    a = list_size // n
    b = list_size % n
    return [list[i*a + (i if i < b else b):(i+1)*a + (i+1 if i < b else b)] for i in range(n)]


class MakeCVDataset:
    def __init__(self, genelabellist, inputfam, validation):
        self.genenum = len(genelabellist)
        self.genelabelindex = []
        self.usingfamlist = []
        for i in inputfam:
            labelindex = np.where(genelabellist==i)[0]
            # Families whose number of genes is fewer than the number of CV aren't used.
            if len(labelindex) >= validation:
                np.random.seed(randomseed)
                # 分割をランダムに行う場合は下のコメントアウトを外す
                #labelindex = np.random.permutation(labelindex)
                self.genelabelindex.append(split_list(list(labelindex), validation))
                self.usingfamlist.append(i)
        #print(self.usingfamlist)
        #print(self.genelabelindex)

    def loaddataset(self, datapath, labelpath):
        print("Loading data...")
        start = time.time()

        self.data = np.load(datapath).astype(np.float32)
        self.label = np.load(labelpath).astype(np.int32)

        print("Complete loading Files.")
        elapsed_time = time.time() - start
        print ("Loding time:{0}".format(elapsed_time) + "[sec]")
        print("Data size :",len(self.data),"*",len(self.data[0]),"*",len(self.data[0][0]))
        print("")

        self.genelen = len(self.data[0])
        self.width = len(self.data[0][0])

    def makeindex(self, val):
        index = np.empty((0), dtype=np.int32)

        ## index
        # train -> 1
        # test -> 2

        trainindex = []
        testindex = []
        for i in range(len(self.usingfamlist)):
            trainindex.extend(self.genelabelindex[i][:val])
            trainindex.extend(self.genelabelindex[i][val+1:])
            testindex.extend(self.genelabelindex[i][val])
        trainindex = [flatten for inner in trainindex for flatten in inner]
        #print(trainindex)
        #print(testindex)

        for i in range(0,self.genenum-1):
            miniindex = np.zeros(self.genenum-i-1).astype(np.int32)
            if i in testindex:                      # testのindex書き込み
                for j in range(i+1, self.genenum):
                    if j in testindex:
                        miniindex[j-i-1] = 2
            elif i in trainindex:                   # trainのindex書き込み
                for j in range(i+1, self.genenum):
                    if j in trainindex:
                        miniindex[j-i-1] = 1
            index = np.concatenate([index, miniindex], axis=0)
            #print("   "*i,list(miniindex))

        return index

    #=========================================================================================

    def splitdataset(self, index, point):     # point = 1:train, 2:test
        partdata = self.data[index==point]
        partlabel = self.label[index==point]

        return partdata, partlabel

    def makeCVdataset(self, index):
        print("Making dataset...")
        start = time.time()
        
        ## Split dataset
        self.train_data = self.data[index==1] 
        self.train_label = self.label[index==1]
        #print(self.train_label.shape)
        self.test_data = self.data[index==2]
        self.test_label = self.label[index==2]
        #print(self.test_label.shape)

        #del self.data
        #del self.label
        #gc.collect()

        ## Increase dataset
        x = np.concatenate([np.zeros(int(self.width/2)),np.ones(int(self.width/2))])
         # Increase train daaset
        self.train_data = np.concatenate([self.train_data,
                                          np.c_[self.train_data[:,:,x==1], self.train_data[:,:,x==0]]], axis=0)
        self.train_label = np.concatenate([self.train_label, self.train_label])
        y = np.random.permutation(len(self.train_data))
        self.train_data = self.train_data[y]
        self.train_label = self.train_label[y]

         # Increase test dataset
        self.test_data = np.concatenate([self.test_data,
                                         np.c_[self.test_data[:,:,x==1], self.test_data[:,:,x==0]]], axis=0)
        self.test_label = np.concatenate([self.test_label, self.test_label])
        y = np.random.permutation(len(self.test_data))
        self.test_data = self.test_data[y]
        self.test_label = self.test_label[y]

        elapsed_time = time.time() - start
        print("Complete making dataset.")
        print("Making time:{0}".format(elapsed_time) + "[sec]")
        print("")

        self.train_data = self.train_data.reshape(len(self.train_data), 1, self.genelen, self.width)
        print("Input train data size :", self.train_data.shape)

        self.test_data = self.test_data.reshape(len(self.test_data), 1, self.genelen, self.width)
        print("Input test data size :", self.test_data.shape)
        print("")

        train = tuple_dataset.TupleDataset(self.train_data, self.train_label)
        test = tuple_dataset.TupleDataset(self.test_data, self.test_label)

        del self.train_data
        del self.test_data
        gc.collect()

        return train, test, self.width

    def getTestdataset(self, index):
        print("Makeing dataset...")
        start = time.time()

        ## Split dataset
        test_data = self.data[index==2]
        test_label = self.label[index==2]
        print(test_label.shape)

        elapsed_time = time.time() - start
        print("Complete make dataset.")
        print("Making time:{0}".format(elapsed_time) + "[sec]")
        print("")

        test_data = test_data.reshape(len(test_data), 1, self.genelen, self.width)
        print("Input test data size :", test_data.shape)
        print("")

        return test_data, test_label, self.width


class MakeUnknownFamiliesDataset:
    def __init__(self, genelabellist, inputfam, testlist):
        self.genenum = len(genelabellist)
        self.genelabellist = genelabellist
        self.testlist = testlist
        self.trainlist = list(set(inputfam)-set(testlist))
        #print(self.testlist)
        #print(self.trainlist)

    def loaddataset(self, datapath, labelpath):
        print("Loading data...")
        start = time.time()

        self.data = np.load(datapath).astype(np.float32)
        self.label = np.load(labelpath).astype(np.int32)

        print("Complete loading Files.")
        elapsed_time = time.time() - start
        print ("Loding time:{0}".format(elapsed_time) + "[sec]")
        print("Data size :",len(self.data),"*",len(self.data[0]),"*",len(self.data[0][0]))
        print("")

        self.genelen = len(self.data[0])
        self.width = len(self.data[0][0])

    def makeindex(self):
        index = np.empty((0), dtype=np.int32)

        ## index
        # train -> 1
        # test -> 2

        trainindex = np.empty((0), dtype=np.int32)
        testindex = np.empty((0), dtype=np.int32)
        for i in self.trainlist:
            trainindex = np.append(trainindex, np.where(self.genelabellist==i)[0])
        for i in self.testlist:
            testindex = np.append(testindex, np.where(self.genelabellist==i)[0])
        trainindex = np.sort(trainindex).tolist()
        testindex = np.sort(testindex).tolist() 
        #print(trainindex)
        #print(testindex)

        for i in range(0,self.genenum-1):
            miniindex = np.zeros(self.genenum-i-1).astype(np.int32)
            if i in testindex:                      # testのindex書き込み
                for j in range(i+1, self.genenum):
                    if j in testindex:
                        miniindex[j-i-1] = 2
            elif i in trainindex:                   # trainのindex書き込み
                for j in range(i+1, self.genenum):
                    if j in trainindex:
                        miniindex[j-i-1] = 1
            index = np.concatenate([index, miniindex], axis=0)
            #print("   "*i,list(miniindex))

        return index

    #=========================================================================================

    def splitdataset(self, index, point):     # point = 1:train, 2:test
        partdata = self.data[index==point]
        partlabel = self.label[index==point]

        return partdata, partlabel

    def makeUFdataset(self, index):
        print("Making dataset...")
        start = time.time()

        ## Split dataset
        self.train_data = self.data[index==1]
        self.train_label = self.label[index==1]
        #print(self.train_label.shape)
        self.test_data = self.data[index==2]
        self.test_label = self.label[index==2]
        #print(self.test_label.shape)

        ## Increase dataset
        x = np.concatenate([np.zeros(int(self.width/2)),np.ones(int(self.width/2))])
         # Increase train dataset
        self.train_data = np.concatenate([self.train_data,
                                          np.c_[self.train_data[:,:,x==1], self.train_data[:,:,x==0]]], axis=0)
        self.train_label = np.concatenate([self.train_label, self.train_label])
        y = np.random.permutation(len(self.train_data))
        self.train_data = self.train_data[y]
        self.train_label = self.train_label[y]

         # Increase test dataset
        self.test_data = np.concatenate([self.test_data,
                                         np.c_[self.test_data[:,:,x==1], self.test_data[:,:,x==0]]], axis=0)
        self.test_label = np.concatenate([self.test_label, self.test_label])
        y = np.random.permutation(len(self.test_data))
        self.test_data = self.test_data[y]
        self.test_label = self.test_label[y]

        elapsed_time = time.time() - start
        print("Complete making dataset.")
        print("Making time:{0}".format(elapsed_time) + "[sec]")
        print("")

        self.train_data = self.train_data.reshape(len(self.train_data), 1, self.genelen, self.width)
        print("Input train data size :", self.train_data.shape)

        self.test_data = self.test_data.reshape(len(self.test_data), 1, self.genelen, self.width)
        print("Input test data size :", self.test_data.shape)
        print("")

        train = tuple_dataset.TupleDataset(self.train_data, self.train_label)
        test = tuple_dataset.TupleDataset(self.test_data, self.test_label)

        del self.train_data
        del self.test_data
        gc.collect()

        return train, test, self.width

    def getTestdataset(self, index):
        print("Makeing dataset...")
        start = time.time()

        ## Split dataset
        test_data = self.data[index==2]
        test_label = self.label[index==2]
        print(test_label.shape)

        elapsed_time = time.time() - start
        print("Complete make dataset.")
        print("Making time:{0}".format(elapsed_time) + "[sec]")
        print("")

        test_data = test_data.reshape(len(test_data), 1, self.genelen, self.width)
        print("Input test data size :", test_data.shape)
        print("")

        return test_data, test_label, self.width
