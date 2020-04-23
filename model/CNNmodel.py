import argparse

import numpy as np
import chainer
from chainer import report, training, Chain, datasets, iterators, optimizers, cuda, Reporter, report_scope
import chainer.functions as F
import chainer.links as L
from chainer.training import extensions
from chainer.datasets import tuple_dataset

import matplotlib.pyplot as plt
plt.switch_backend('agg')


from chainer import serializers
from chainer.cuda import to_gpu
from chainer.cuda import to_cpu


class TestModeEvaluator(extensions.Evaluator):
    
    def evaluate(self):
        model = self.get_target('main')
        model.train = False
        ret = super(TestModeEvaluator, self).evaluate()
        model.train = True
        return ret

class CNN(Chain):
    def __init__(self, output_ch1,output_ch2, filter_height, filter_width, n_units, n_label):
        super(CNN, self).__init__(
            conv1 = L.Convolution2D(1, output_ch1, (filter_height, filter_width)),
            bn1 = L.BatchNormalization(output_ch1),
            conv2 = L.Convolution2D(output_ch1, output_ch2, (filter_height, 1)),
            bn2 = L.BatchNormalization(output_ch2),
            fc1 = L.Linear(None, n_units),
            fc2 = L.Linear(None, n_label))

    def __call__(self, x, train=True):
        h1 = F.max_pooling_2d(F.relu(self.bn1(self.conv1(x))), ksize=(10,1), stride=8)
        h2 = F.max_pooling_2d(F.relu(self.bn2(self.conv2(h1))), ksize=(14,1), stride=8)
        h3 = F.dropout(F.relu(self.fc1(h2)), ratio=0.5)
        y = self.fc2(h3)
        return y
        
class MyClassifier(Chain):
    def __init__(self, predictor, lastlayer):
        super(MyClassifier, self).__init__(predictor=predictor)
        self.lastlayer = lastlayer

    def __call__(self, x, t):
        y = self.predictor(x)
        if self.lastlayer == 1:     # The number of last layer units = 1
            loss = F.sigmoid_cross_entropy(y, t.reshape(len(t),1))
            accuracy = F.binary_accuracy(y, t.reshape(len(t),1))
        else:                       # The number of last layer units = 2
            loss = F.softmax_cross_entropy(y, t)
            accuracy = F.accuracy(y, t)
        summary = F.classification_summary(y, t, beta = 1.0)
        precision = summary[0]
        recall = summary[1]
        f_value = summary[2]
        reporter = Reporter()
        observer = object()
        reporter.add_observer('f_value:', observer)
        observation={}    
        with reporter.scope(observation):
            reporter.report({'x':f_value},observer)    
        report({'loss': loss,
                'accuracy': accuracy, 
                'precision': precision,
                'recall': recall,
                'f_value': f_value}, self)
        report(dict(('f_value_%d' % i, val) for i, val in enumerate(f_value)), self)
        #report(dict(('precision_%d' % i, val) for i, val in enumerate(precision)), self)
        #report(dict(('recall_%d' % i, val) for i, val in enumerate(recall)), self)
        return loss 


class Runchainer:
    def __init__(self, batchsize, outdir, output_ch1, output_ch2, filter_height, width, n_units, n_label):
        self.batchsize = batchsize
        self.outdir = outdir
        self.output_ch1 = output_ch1
        self.output_ch2 = output_ch2
        self.filter_height = filter_height
        self.width = width
        self.n_units = n_units
        self.n_label = n_label

    def learn(self, train, test, gpu, epoch):

        self.epoch = epoch

        ## Set up model =============================================
        model = MyClassifier(CNN(self.output_ch1, self.output_ch2, self.filter_height, \
                                 self.width, self.n_units, self.n_label), self.n_label)
        # Set up GPU ------------------------------------------
        if gpu >= 0:
            chainer.cuda.get_device(gpu).use()
            model.to_gpu()

        ## Set up Optimizer =========================================
        optimizer = optimizers.Adam()
        optimizer.setup(model)

        ## Set up iterator ==========================================
        train_iter = iterators.SerialIterator(train, self.batchsize)
        test_iter = chainer.iterators.SerialIterator(test, self.batchsize,
                                                 repeat=False, shuffle=False)

        ## Updater ==================================================
        updater = training.StandardUpdater(train_iter, optimizer, device=gpu)

        ## Trainer ==================================================
        trainer = training.Trainer(updater, (self.epoch, 'epoch'), out=self.outdir)

        ## Evaluator ================================================
        trainer.extend(TestModeEvaluator(test_iter, model, device=gpu))
        #trainer.extend(extensions.ExponentialShift('lr', 0.5),
        #               trigger=(25, 'epoch'))
        trainer.extend(extensions.dump_graph('main/loss'))
        trainer.extend(extensions.snapshot(filename='snapshot'))
        trainer.extend(extensions.snapshot_object(model.predictor, filename='model_epoch-{.updater.epoch}'))
        trainer.extend(extensions.LogReport(log_name='log'))
        trainer.extend(extensions.PrintReport(
            ['epoch', 'main/accuracy', 'validation/main/accuracy', \
             'main/f_value_0', 'validation/main/f_value_0', 'elapsed_time']))

        trainer.extend(extensions.PlotReport(
            ['main/loss', 'validation/main/loss'],
            x_key='epoch', file_name='loss.png'))
        trainer.extend(extensions.PlotReport(
            ['main/accuracy', 'validation/main/accuracy'],
            x_key='epoch', file_name='accuracy.png'))

        trainer.extend(extensions.ProgressBar())

        trainer.run()

    def predict(self, test_data, test_label, predictorpath):
        
        predictor = CNN(self.output_ch1, self.output_ch2, self.filter_height, \
                        self.width, self.n_units, self.n_label)
        serializers.load_npz(predictorpath, predictor)

        TP = 0
        FN = 0
        FP = 0
        TN = 0

        for i in range(len(test_label)):
            predictor.to_cpu()
            #y = model(test_data[i].reshape(1,1,1200,16)).data.argmax(axis=1)#[0]
            ## 最終層を1層にした場合は出力が単一で実数となる.
            ## そのため、出力値が負の場合は0, 正の場合は1になる.
            y = predictor(test_data[i].reshape(1,1,-1,self.width)).data.argmax(axis=1)[0]

            if test_label[i] == 0:
                if y == 0:
                    TP += 1
                else:
                    FN += 1
            else:
                if y == 0:
                    FP += 1
                else:
                    TN += 1

        try:
            precision = TP/(TP+FP)
        except ZeroDivisionError:
            precision = 0
        try:
            recall = TP/(TP+FN)
        except ZeroDivisionError:
            recall = 0
        try:
            Fvalue = 2*precision*recall/(precision+recall)
        except ZeroDivisionError:
            Fvalue = 0
        try:
            Accuracy = (TP+TN)/(TP+FN+FP+TN)
        except ZeroDivisionError:
            Accuracy = 0
        
        text = "True Positive : " + str(TP) + "\n" + \
               "False Negative : " + str(FN) + "\n" + \
               "False Positive : " + str(FP) + "\n" + \
               "True Negative : " + str(TN) + "\n\n" + \
               "Accuracy : " + str(Accuracy) + "\n" + \
               "Precision : " + str(precision) + "\n" + \
               "Recall : " + str(recall) + "\n" + \
               "Fvalue : " + str(Fvalue) + "\n"

        return text
