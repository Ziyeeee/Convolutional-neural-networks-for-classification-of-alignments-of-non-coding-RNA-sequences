# -*- coding: utf-8 -*-

import argparse
import sys
import subprocess
import os.path

import numpy as np
import chainer
from chainer import report, training, Chain, datasets, iterators, optimizers, cuda, Reporter, report_scope
import chainer.functions as F
import chainer.links as L
from chainer.training import extensions
from chainer.datasets import tuple_dataset

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import time
import gc

sys.path.append("./model")
import DatasetMaker as Dataset
import CNNmodel as Model

from chainer import serializers

class LearningMode:
    def __init__(self, args):
        self.dataset = args.dataset
        self.label = args.label
        self.batchsize = args.batchsize
        self.epoch = args.epoch
        self.validation = args.validation
        self.gpu = args.gpu
        self.out = args.out
        self.resume = args.resume

        ## Process genelabel data
        genelabellist = []
        for line in open(args.genelabel, 'r'):
            lines = line[:-1].split(":")    # Separater is ":".
            genelabellist.append(int(lines[1]))
        self.genelabellist = np.asarray(genelabellist)

        self.inputfam = np.unique(self.genelabellist)

        ##### hyperparameter ###################################################
        self.output_ch1 = 64        # number of 1st convolutional layer's filters
        self.output_ch2 = 128       # number of 2nd convolutional layer's filters
        self.filter_height = 15     # filter size of convolutional layer
        self.n_units = 1632         # number of middel layer's units
        self.n_label = 2            # number of last layer's units
        ########################################################################

        ## For unknown families mode
        self.testlist = [0,1,5]     # ID of families that be used in test
                                    # 0 : snRNA
                                    # 1 : snoRNA C/D-box
                                    # 2 : snoRNA H/ACA-box
                                    # 3 : miRNA
                                    # 4 : YRNA
                                    # 5 : tRNA
                                    # (6 : SCARNA, 7 : Vault RNA, 8 : 5S rRNA, 9 : others)
        genelist = ["snRNA", "snoRNA_C/D", "snoRNA_H/ACA", "miRNA", "YRNA", "tRNA",
                    "SCARNA", "Vault_RNA", "5S_rRNA", "others"]
        self.inputfamname = []           # list of inputted gene families' name
        for i in self.inputfam:
            self.inputfamname.append(genelist[i])
        self.traingenelist = []
        self.testgenelist = []
        for i in self.inputfam:
            if i in self.testlist:
                self.testgenelist.append(genelist[i])
            else:
                self.traingenelist.append(genelist[i])

    def output_parameterlog(self):
        text = "*Hyperparameter\n" + \
               " 1st conv filter num : " + str(self.output_ch1) + "\n" + \
               " 2nd conv filter num : " + str(self.output_ch2) + "\n" + \
               " Conv filter height : " + str(self.filter_height) + " * [feature vector size]\n" + \
               " Middle layer units num : " + str(self.n_units) + "\n" + \
               " Last layer units num : " + str(self.n_label) + "\n" + \
               " Minibatch size : " + str(self.batchsize) + "\n" + \
               " epoch : " + str(self.epoch) + "\n" + \
               "\n" + \
               "*GPU id : " + str(self.gpu) + "\n" + \
               "\n" + \
               "*Dataset\n" + \
               " input data : " + self.dataset + "\n" + \
               " input label : " + self.label + "\n" + \
               " used gene families : " + ', '.join(self.inputfamname) + "\n"
        ## For unknown families mode
        if self.validation <= 0:
            text += " -train gene family : " + ', '.join(self.traingenelist) + "\n" + \
                    " -test gene family  : " + ', '.join(self.testgenelist) + "\n"
                   
        if os.path.isdir(self.out)==False:
            subprocess.call(["mkdir", self.out])
        f = open(self.out+'/parameterlog', 'w')
        f.write(text)
        f.close()

    def cross_validation(self, validation):
        CV = Dataset.MakeCVDataset(self.genelabellist, self.inputfam, self.validation)
        CV.loaddataset(self.dataset, self.label)

        val = 9
        valend = validation
        while val < valend:
            print("Cross Validation No.", val, "=====================================")
            index = CV.makeindex(val)
            train, test, width = CV.makeCVdataset(index)

            outdir = self.out + "/validation" + str(val)
            RC = Model.Runchainer(self.batchsize, outdir, self.output_ch1, self.output_ch2,
                                  self.filter_height, width, self.n_units, self.n_label)
            RC.learn(train, test, self.gpu, self.epoch)

            del RC
            gc.collect()
            print("")

            val += 1

    def unknown_families(self):
        print("Train genes :", self.traingenelist)
        print("Test genes  :", self.testgenelist)
        print("")

        UF = Dataset.MakeUnknownFamiliesDataset(self.genelabellist, self.inputfam, self.testlist)
        UF.loaddataset(self.dataset, self.label)

        index = UF.makeindex()
        train, test, width = UF.makeUFdataset(index)

        outdir = self.out
        RC = Model.Runchainer(self.batchsize, outdir, self.output_ch1, self.output_ch2,
                              self.filter_height, width, self.n_units, self.n_label)
        RC.learn(train, test, self.gpu, self.epoch)


class PredictingMode(LearningMode):
    def __init__(self, args):
        super().__init__(args)
        self.args = args

    def predictor(self):
        if self.args.validation > 0:
            val = 0

            CV = Dataset.MakeCVDataset(self.genelabellist, self.inputfam, self.validation)
            CV.loaddataset(self.dataset, self.label)

            index = CV.makeindex(val)
            test_data, test_label, width = CV.getTestdataset(index)

        else:
            UF = Dataset.MakeUnknownFamiliesDataset(self.genelabellist, self.inputfam, self.testlist)
            UF.loaddataset(self.dataset, self.label)

            index = UF.makeindex()
            test_data, test_label, width = UF.getTestdataset(index)
            
        outdir = self.out
        RC = Model.Runchainer(self.batchsize, outdir, self.output_ch1, self.output_ch2,
                              self.filter_height, width, self.n_units, self.n_label)
        predictresult = RC.predict(test_data, test_label, self.args.predictor)

        f = open(self.args.out+'/predict_result', 'w')
        f.write(predictresult)
        f.close()

def main():
    parser = argparse.ArgumentParser(description='Chainer ncRNA classification:')
    parser.add_argument('--dataset', '-d', default="./testdata/ncRNApairdataDAFS_test.npy",
                        help='The dataset to use: numpy file (hoge.npy) of aligned 2 ncRNAs pair data.')
    parser.add_argument('--label', '-l', default="./testdata/ncRNApairlabel_test.npy",
                        help='The label to use: numpy file (huga.npy) corresponding to the dataset.')
    parser.add_argument('--genelabel', '-gl', default="./testdata/genelabel_6families.txt",
                        help='The gene annotation file (colon separated).')
    parser.add_argument('--batchsize', '-b', type=int, default=128,
                        help='Number of images (aligned 2 ncRNAs pair data) in each mini-batch.')
    parser.add_argument('--epoch', '-e', type=int, default=25,
                        help='Number of sweeps over the dataset to train.')
    parser.add_argument('--validation', '-v', type=int, default=10,
                        help='Number of cross validation (if v<=0, run in the unknown families mode).')
    parser.add_argument('--gpu', '-g', type=int, default=-1,
                        help='GPU ID (negative value indicates CPU)')
    parser.add_argument('--out', '-o', default='result',
                        help='Directory to output the result')
    parser.add_argument('--resume', '-r', default='',
                        help='Resume the training from snapshot')
    parser.add_argument('--predictor', '-p', default='',
                        help='Learned model file to predict data')
    args = parser.parse_args()

    if args.predictor:
        print("### Predicting mode ###\n")
        PM = PredictingMode(args)
        PM.predictor()

    else:
        print("### Learning mode ###\n")
        print('GPU: {}'.format(args.gpu))
        print('# Minibatch-size: {}'.format(args.batchsize))
        print('# epoch: {}'.format(args.epoch))
        print('# validation num: {}'.format(args.validation))
        print('')

        LM = LearningMode(args)
        LM.output_parameterlog()
        
        if args.validation > 0:
            print(" ## Cross-validation mode ##\n")
            LM.cross_validation(args.validation)

        else:
            print(" ## Unknown families mode ##\n")
            LM.unknown_families()

if __name__ == '__main__':
    main()
