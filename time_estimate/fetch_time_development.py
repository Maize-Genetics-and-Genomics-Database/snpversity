import csv
import statsmodels.formula.api as sm
import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np
from sklearn import neighbors
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import AdaBoostRegressor
import os.path
from sklearn.externals import joblib
import math
from mpl_toolkits.mplot3d import Axes3D
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
import random


n_neighbors = 3
KNN_WEIGHTS = 'uniform'
PLOT_STOCKS_POINTS = 100.0
PLOT_RANGES_POINTS = 10.0
MODEL_TYPE = 'regression'  # 'regression' OR 'ada' OR 'knn'
ADA_N_ESTIMATORS = 200
VALIDATION_TYPE = 'mae'
DEV_IDENTIFIER = '_development'

def parseCSV(name):
    xs = []
    ys = []
    zs = []
    with open(name) as allzeagbs:
        reader = csv.reader(allzeagbs, delimiter=',')
        for row in reader:
            x = int(row[0])
            y = int(row[1])
            z = int(round(float(row[2])))
            xs.append(x)
            ys.append(y)
            zs.append(z)
            co_ord = (x,y,z)
            #print('(x,y,z) = {0},{1},{2}'.format(x,y,z))
    return xs,ys,zs


def plot(dataset, model_type='knn', point_selection='random'):
    taxa,ranges,times = parseCSV(dataset+DEV_IDENTIFIER+".csv")
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('# Stocks')
    ax.set_ylabel('Ranges (bp)')
    ax.set_zlabel('Time (s)')
    if point_selection == 'step':
        points_generated = [(s, r) for s in xrange(min(taxa), max(taxa), int(math.ceil((max(taxa) - min(taxa)) / PLOT_STOCKS_POINTS)))
                            for r in xrange(min(ranges), max(ranges), int(math.ceil((max(ranges) - min(ranges)) / PLOT_RANGES_POINTS)))]
    else:
        points_generated = [(random.randrange(min(taxa), max(taxa)), random.randrange(min(ranges), max(ranges)))
                            for i in xrange(0,int(PLOT_STOCKS_POINTS * PLOT_RANGES_POINTS))]

    if model_type in ['knn', 'ada']:
        model = load(dataset, model_type)
        time_predictions = [predict(model,s,r, model_type) for s,r in points_generated]
        plt.title(dataset + " ({0})".format(model_type))
    else:
        model = load(dataset)
        time_predictions = [predict(model,s,r, 'regression') for s,r in points_generated]
        plt.title(dataset + " (multivariable regression)")

    stocks_generated = [s for s,r in points_generated]
    ranges_generated = [r for s,r in points_generated]
    ax.scatter(taxa,ranges,times,c='r', marker='o', label='Actual')
    ax.scatter(stocks_generated,ranges_generated,time_predictions,c='b',marker='^', label='Predicted')
    plt.legend(loc='upper left')
    plt.show()


def construct(dataset,stocks=[],ranges=[],times=[]):
    if not stocks or not ranges or not times:
        stocks,ranges,times = parseCSV(dataset+DEV_IDENTIFIER+".csv")
    df = pd.DataFrame({"t":times, "s":stocks, "r":ranges})
    return sm.ols("t ~ s + r", data=df).fit() # Time ~ Stocks, Range


def predict(model, s, r, model_type='regression'):
    if model_type in ['knn', 'ada']:
        t_ = model.predict(np.asarray([s,r]).reshape(1, -1))
        #print(t_[0])
        return abs(t_[0])
    else:
        t = model.params.Intercept + (s*model.params.s) + (r*model.params.r)
        #print(t)
    return abs(t)


def scikit_construct(dataset, mode='knn', X=np.array([]), y=np.array([])):
    if not X.any() or not y.any():
        taxa,ranges,times = parseCSV(dataset+DEV_IDENTIFIER+".csv")
        X = np.asarray(zip(taxa,ranges))
        y = np.asarray(times)
    if mode == 'knn':
        knn = neighbors.KNeighborsRegressor(n_neighbors, weights=KNN_WEIGHTS)
        model = knn.fit(X, y)
    elif mode == 'ada':
        rng = np.random.RandomState(1)
        ada = AdaBoostRegressor(DecisionTreeRegressor(max_depth=4), n_estimators=ADA_N_ESTIMATORS, random_state=rng)
        model = ada.fit(X, y)
    return model


def load(dataset, mode_type='regression'):
    fname = dataset + '_development.pickle'
    if mode_type not in ['knn', 'ada']:
        from statsmodels.iolib.smpickle import load_pickle
        if os.path.isfile(fname):
            return load_pickle(fname)
        else:
            print("Pickled file " + fname + " not found!")
            print("Constructing model...")
            model = construct(dataset)
            print("Serializing...")
            model.save(fname)
            return model
    else:
        fname = dataset + '_{0}.pickle'.format(mode_type)
    if os.path.isfile(fname):
        return joblib.load(fname)
    else:
        print("Pickled file " + fname + " not found!")
        print("Constructing model...")
        model = scikit_construct(dataset, mode_type)
        print("Serializing...")
        joblib.dump(model, fname)
        return model


def load_train(dataset,stocks,ranges,times, mode_type='regression'):
    if mode_type not in ['knn', 'ada']:
        model = construct(dataset,stocks,ranges,times)
    else:
        X = np.asarray(zip(stocks,ranges))
        y = np.asarray(times)
        model = scikit_construct(dataset,mode_type,X,y)
    return model


def model_validation(dataset, model_type='regression', validation_type='mse', verbose=False):
    """
    Calculates mean error of model's leave-one-out
    :param dataset: File path of CSV file containing data points.
    :param model_type: regression, knn, ada
    :param validation_type: mse=Mean Squared Error, mae=Mean Absolute Error
    :param verbose: Print out actual and predicted results
    :return:
    """
    loo = LeaveOneOut()
    taxa,ranges,times = parseCSV(dataset+".csv")
    # For experimentation between test/train taxa_dev,ranges_dev,times_dev = parseCSV(dataset+DEV_IDENTIFIER+".csv")
    print "Evaluation of {0} using {1} and {2} real datapoints".format(dataset,model_type,str(len(taxa)))
    if verbose:
        print "# Stocks, Range (bp), Actual time, Predicted time"
    X = np.asarray(zip(taxa,ranges))
    y = np.asarray(times)
    time_true = []
    time_pred = []
    for train,test in loo.split(X):
        # Generate train dataset and leave one out for testing
        taxa_train = [taxa[i] for i in train]
        ranges_train = [ranges[i] for i in train]
        times_train = [times[i] for i in train]
        model = load_train(dataset,taxa_train,ranges_train,times_train,model_type)
        s_test = taxa[test] # s = stocks = taxa
        r_test = ranges[test] # r = range
        time_true.append(times[test]) # true = actual
        pred = predict(model,s_test,r_test,model_type)
        time_pred.append(pred)
        if verbose:
            print '{0},{1},{2},{3}'.format(s_test,r_test,times[test],pred)
    if validation_type == 'mae':
        error = mean_absolute_error(time_true,time_pred)
    else:
        error = mean_squared_error(time_true, time_pred)
    print ('Mean of {1}: {0}'.format(mean(time_true),dataset))
    print ('Standard Deviation of {1}: {0}'.format(np.std(time_true),dataset))
    return error


def mean(numbers):
    return float(sum(numbers) / max(len(numbers),1))



def main():
    global MODEL_TYPE
    if len(sys.argv) == 4:
        # Predict
        dataset = sys.argv[1]
        s = int(sys.argv[2])
        r = int(sys.argv[3])
        if dataset == 'ZeaHM321_raw':
            MODEL_TYPE = 'knn'
        model = load(dataset, MODEL_TYPE)
        t = predict(model,s,r, MODEL_TYPE)
        print t
    elif len(sys.argv) == 2:
        dataset = sys.argv[1]
        if dataset == 'all':
            for d in ['AllZeaGBSv27public20140528','ZeaGBSv27publicImputed20150114','ZeaHM321_LinkImpute','ZeaHM321_raw']:
                mse = model_validation(d,MODEL_TYPE,VALIDATION_TYPE)
                print '{2} of {0}: {1}'.format(d,mse,VALIDATION_TYPE)
                plot(d, MODEL_TYPE,'uniform')
                #plot(d, 'knn','uniform')
                #plot(d, 'ada')
        else:
            mse = model_validation(dataset,MODEL_TYPE,VALIDATION_TYPE,True)
            print '---------SUMMARY---------'
            print '{2} of {0} using {3} model: {1}'.format(dataset,mse,VALIDATION_TYPE,MODEL_TYPE)
            plot(dataset, MODEL_TYPE,'step')
    else:
        sys.exit("Please specify dataset, Stock-Count, Range of positions. Example:\n"
                 "`python2.7 fetch_time.py AllZeaGBSv27public20140528 100 128121`")


if __name__ == "__main__":
    main()