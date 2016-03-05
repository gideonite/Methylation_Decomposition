# import scipy.io
import pdb
from matplotlib.pylab import *
import numpy as np

def main():

    dists = []
    bins = [0,2,4,6,8,12,16,24,32,48,64,128,256,512,1024,10000]
    bins = np.asarray(bins)
    totals = [np.zeros(2) for i in range(len(bins)-1)]
    cts = [np.zeros((2,2)) for i in range(len(bins)-1)]

    print bins
    for cell in [1,2,3]:
        line_prev = ['init']
        filename = '../data/all_truth/cell%d.out' %(cell,)
        with open(filename) as f:
            for line in f:
                line_this = line.split('\t')
                assert(len(line_this) == 6)
                if line_this[0] == line_prev[0]:
                    # assume monotonic position (in either direction)
                    dist = abs(int(line_prev[3])- int(line_this[3]))
                    prev = int(float(line_prev[4]))
                    this = int(float(line_this[4]))
                    dists.append(dist)
                    idx = np.argmax(bins>=dist)-1
                    totals[idx][prev]+=1
                    cts[idx][prev][this]+=1
                line_prev = line_this
    freqs = [ct/np.transpose(np.vstack([tt,tt])) for ct,tt in zip(cts,totals)]
    scipy.io.savemat('../results/count_bigrams.mat', {'bins':bins,'freqs':np.vstack(freqs)})
    pdb.set_trace()
    figure()
    events, edges, patches = hist(dists, bins=bins)
    xlim((0,2000))
    print events
    show()




if __name__ == '__main__':
    main()