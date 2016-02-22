#!/usr/bin/env python
import sys
import numpy as np
import scipy.io
import argparse
import os
import re
import os.path
from operator import itemgetter

'''
    Sets up command line info
'''
def parse_args(argv):
    parser = argparse.ArgumentParser(prog="doCollection",
                                     description=('Collection over all graphs'))
    parser.add_argument('directory',help='directory for the mat files')
    parser.add_argument('output',help='output for the new mat file')

    args= parser.parse_args(argv[1:])
    return args.directory, args.output


'''
    This will go through a directory and get the Phase lag and Pop Correl
    of each file. I am assuming the directory you are looking through
    will only contain postproccesed data and the values will be in the file
'''
def collect_plpc(dir_path):
    data = []
    for f in os.listdir(dir_path):
        ##Only care about mat files 
        if f.endswith('.mat'):
            file_loc = os.path.join(dir_path,f)
            values = scipy.io.loadmat(file_loc)
            phase_lag = values['phase_lag']
            pop_correlation = values['pop_correlation']

            ##This gets the degree so we can arrange data since f is not in order
            string_pieces = re.split(r'[_]+|idegree',f)
            degree = string_pieces[3]

            ##Sends data as tuple for easier graphing
            data.append((degree,phase_lag,pop_correlation))

    ##Now we sort the data according to the degree to make it ordered
    ordered_data =  sorted(data,key=itemgetter(0))

    ##Change list to np.array so that it can be saved as .mat
    data = np.asarray(ordered_data,dtype=object)
    return data

def main(argv=None):
    if argv is None:
        argv = sys.argv
    (directory, out_fn) = parse_args(argv)

    phase_lag_pop_corr = collect_plpc(directory)

    scipy.io.savemat(out_fn,
                     mdict = {'phase_lag_pop_corr':phase_lag_pop_corr}, oned_as = 'column')
if __name__ == '__main__':
    status = main()
    sys.exit(status)
    
            
