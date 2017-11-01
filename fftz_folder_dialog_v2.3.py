#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python 2.7
# TODO make description read good
# TODO x and y signals
# TODO time domain signal output solution
""" FFT and PSD calculator for NanoTracker optical trap data
    Read NT force-save.txt files
    Take fft and calc psd
    Write time domain and PSD spectrums to csv files

0.2.1:
    - Filter problems fix:
        - rm the bandpass_ifft filter that hates >25KHz signal
        - implemented a order 3 butterworth highpass filter (pass = >1Hz)
    - experimental gui folder dialog to select root data folder

0.2.2
    - bugfix, outputfile path wrong on Windows due to dirsep char from 
      gui folder dialog
"""
### Imports and defs {{{
# -------------------------------------------------
import argparse
import csv
import numpy as np
from glob import glob
import os
import sys
import pandas as pd
import scipy.signal
class CommentedFile:
    """ Skips all lines that start with commentstring
    """
    def __init__(self, f, commentstring="#"):
        self.f = f
        self.commentstring = commentstring
    def next(self):
        line = self.f.next()
        while line.startswith(self.commentstring):
            line = self.f.next()
        return line
    def __iter__(self):
        return self
def dirsep():
    """Check OS and set appropriate dir separator 
    """ 
    if sys.platform.lower().startswith('win'):
        dirsep = '\\'
    else:
        dirsep = '/'
    return dirsep
def get_cwd():
    ### Check OS switch 
    if sys.platform.lower().startswith('win'):
        cdir = os.getcwd().split('\\')[-1] # Windows
    else:
        cdir = os.getcwd().split('/')[-1] # *.nix
    return cdir
def freq_calc(N, fs):    
    """ Where N is length of time-series signal (1.5 million data points for 15 secs RAMP design)
    Get the freq values for psd (assumes all scans are same length)
    Index + and - components of fz, start at 1 to lose DC offset 
    """
    dt = float(1/fs)
    if N % 2 == 0: ## %=0 means it is an even number of data points. #if even, #F Neq is at index N/2 - if odd, #F Neq is at index N/2 + 1
        ind = np.arange(1,N/2)
    else:
        ind = np.arange(1,N/2+1)   ##ind = index
    # Get corresponding fft freqs 
    freq = np.fft.fftfreq(N, dt)
    # drop dc offset and neg freqs
    freq = freq[ind]
    return freq
def read_forcesave(f, col=2):
    """ Reads NanoTracker force-save.txt
    Depends on: class CommentedFile
    In:
        f: force-save file to be opened
        col: column of f to be extracted, 0=x1, 3=x2, 1=Y1, 4=Y2, 2=z1, 5=z2  ##salim
    Returns:
        sig: 1D time series signal
    """
    # Probs a more robust way to do this using with open() as :
    csv_file = csv.reader(CommentedFile(open(f, "rb")),
                          delimiter=' ')
    
    # This loop is to skip the Position Ramp segment, the next segment, 'Pause'
    # is preceeded by a blank line in the force-save-*.txt
    line = csv_file.next()
    while line: # While the next line is not blank 
        line = csv_file.next()
    
    # Pull out the signal into a list of strings
    sig=[]
    for row in csv_file:
        try:
            sig.append(row[col])
        except IndexError:
            print("IndexError. Something is up here: ")
            print(row)
    
    # List of strings >> list of floats
    sig = map(float, sig)   #float = numbers with decimal points
    return sig


# NOTE test: export filtered, unwindowed sig
def psd_powerplay(X, dt):  #X is a varaible name 
    #TODO only need to calc freq once, should be seperate code
    # OR calc each freq, then compare to ensure same, else return warning
    """Does some PSD stuff
    Inputs:
        X: 1D time series, list of floats
        dt: Sampling rate (1/fs) 
    Returns:
        psd: power spectral density of fft(X)
    Relies on: bandpass_ifft
    """
    N = len(X) # Number of samples
    
    # Baseline corrected signal (shifts signal towards y=0)
    X = X - np.mean(X)

    # Highpass filter -- rm signals below cutoff freq
    if args.filter_on == True:
        # Salim wants nowindow filtered sig as an output
        # NOTE Maybe we could window after butterworth? Does it introduce 
        # discontinuity artifacts?
        X_output = butter_highpass_filter(X, highcut=args.cf_low, fs=1/dt, order=3)

        # Windowing for the signal to be used for PSD calc
        X_psd = X * np.hanning(N)
        X_psd = butter_highpass_filter(X, highcut=args.cf_low, fs=1/dt, order=3)
    else:
        X_output = X
        X_psd = X * np.hanning(N)
    
    # fft transform the data
    fx = np.fft.rfft(X_psd)
    # Index + and - components of fz, start at 1 to lose DC offset 
    if N % 2 == 0:
        ind = np.arange(1,N/2)
    else:
        ind = np.arange(1,N/2+1)

    # Get that PSD, son!
    psd = abs(fx[ind])**2    #to calculate the Neq frequency and PSD, fx is fouriur transofromed data and abs = absolute, **= square
    return psd, X_output
def scan_transformation(infile, scan_name):
    """ Run the above functions and store psd for each scan in a dict (psd_d)
        Also store extracted raw signal in a dict
    """
    if len(infile) > 0:   #to check if the directory that it is in, there is a force save file. infile= number of files that are force.save.txt
        colnames = []     #naming
        infile.sort()
        for i in range(len(infile)):   #for= for each force.save file
            colnames.append(''.join((scan_name, '-', str(i+1))))
            sig = read_forcesave(infile[i])
            psd, sig_corrected = psd_powerplay(sig, dt)
            psd_d[colnames[i]] = psd 
            sig_d[colnames[i]] = sig_corrected
            print(' '.join(("Finished:", colnames[i])))
def dict_to_df(d):
    """ Make df from dict with Hz as 1st col and alpha-num order after
        This fun is mostly about formatting colnames
    """
    df = pd.DataFrame.from_dict(d)
    # Case-insensitive alpha sorting
    df = df.reindex_axis(sorted(df.columns, key=lambda x: x.lower()), axis=1)
    # Formatting cols
    cols = df.columns.tolist()
    for i in range(len(cols)):
        cols[i] = cols[i].title()
    df.columns = cols
    if 'Hz' in cols:
        # Move Hz to first col
        cols.insert(0, cols.pop(cols.index('Hz')))
        df = df.reindex(columns=cols)
    return df
def get_match_val(f, match):
    search = open(f, "r")
    for line in search.readlines():
        if match in line.split():
            return line.split()[-1]
            break
def get_params(f):
    """ Returns: t, fs, date (ndarray, float, string)
        Maybe use date for start of file name?
    """
    # TODO num segments -> num loops (only rel for old ramps)
    T = float(get_match_val(f, "settings.segment.1.duration:"))    # t= duration in secs
    N = int(get_match_val(f, "settings.segment.1.num-points:"))     #N= data points
    fs = N/T
    dt = 1/fs
    t = np.arange(dt, T+dt, dt)  #to generate X axis time point, dt= is the time resolution to make it start and end at the same recording times 
    date = get_match_val(f, "bead-id:")
    date = date.split('-')[0] # date is at front of bead-id string
    return t, fs, date
def butter_highpass(highcut, fs, order=3):
    """ Design a digital highpass Butterworth filter, return filter coefficents
    highcut, fs in Hz
    Shouldn't have to convert to rad/s. As long as freqs are expressed with 
    consistent units, the scaling should take care of the normalization
    from scipy cookbook
    """
    nyq = 0.5 * fs
    norm_highcut = highcut / nyq
    b, a = scipy.signal.butter(order, norm_highcut, btype='high')
    return b, a
def butter_highpass_filter(data, highcut, fs, order=3):
    """ Return Butterworth HP filtered time series
    highcut, fs in Hz
    from scipy cookbook
    """
    b, a = butter_highpass(highcut, fs, order=order)
    y = scipy.signal.lfilter(b, a, data)
    return y
def folder_dialog():
    """ Prompt user to select a dir, return its path
    """
    import Tkinter, tkFileDialog
    root = Tkinter.Tk()
    root.withdraw()
    dirpath = tkFileDialog.askdirectory(parent=root,initialdir="./",
        # dirpath will be to dir that user IS IN when they click confirm
        title='Please select your experiment directory (be IN this folder)')
    return dirpath


# NOTE test OK but still new
def walk_get_params(path):
    """ Read scan params from first force-save encountered and break, assumes 
    all files have the same parameters
    This is a quick fix, and probably bad form but should do alright in a pinch
    """
    success = False
    for subdir, dirs, files in os.walk(path):
        if success == False:
            os.chdir(subdir)
            infile = glob('force-save*.txt')
            if len(infile) > 0:
                print('Reading Scan Parameters from dir: %s' %subdir)
                t, fs, date = get_params(infile[0])
                success == True
                break
    # used by psd_powerplay(N, dt)
    dt = float(1/fs)
    freq = freq_calc(len(t), fs)
    return t, dt, freq
#------------------------------------------------- }}}
### SCRIPT INFO
version = '0.2.3'
day = '2017-06-09'
codename = 'SwampMonkey alpha' # xxSnowPrincessxx
print("Version %s (%s) -- \"%s\"" %(version, day, codename) )


### SCRIPT ARGUMENTS (Lowpass filtering below 1Hz) salim
parser = argparse.ArgumentParser()
parser.add_argument('-cl', '--cf_low', type=float, default=1.,
    help='Low cutoff freq for bandpass filter (Hz)')
parser.add_argument('--nofilter', dest='filter_on', action='store_false')
parser.add_argument('--filter', dest='filter_on', action='store_true')
parser.set_defaults(filter_on=True) #DEF: exported sig will also be filtered
# NOTE testing
parser.add_argument('-x', '--xfft', dest='run_x_fft', action='store_true',
    help='Run FFT on X signal (experimental)')
parser.set_defaults(run_x_fft=False)
parser.add_argument('-y', '--yfft', dest='run_y_fft', action='store_true',
    help='Run FFT on Y signal (experimental)')
parser.set_defaults(run_y_fft=False)
# /test
args = parser.parse_args()



### VARIABLE ASSIGNMENT (this is a discription of all variables I am using) salim
# NOTE Apparently when we import through folder dialog sepchar is '/'
path = folder_dialog() # user selects starting folder (QT)
scan_name = path.split('/')[-1]
outfile = ''.join((scan_name, "_psd.csv")) 
outfile = '/'.join((path, outfile))
outfile = os.path.normpath(outfile) # redundant but let's be safe, I think
outsig = ''.join((scan_name, "_sig.csv")) 
outsig = '/'.join((path, outsig))
outsig = os.path.normpath(outsig) # similar reduncancy to above

#find first force save.txt.file, we can find and read and calculate
t, dt, freq = walk_get_params(path)


### Init data dictionaries
psd_d = {}
psd_d['Hz'] = freq
sig_d = {}
sig_d['Time'] = t
    
### Attempted os looping
for subdir, dirs, files in os.walk(path):
    os.chdir(subdir)
    scan_name = get_cwd()
    print(' '.join(("Entering", scan_name)))
    scan_transformation(glob('force-save*.txt'), scan_name)

### Make and export dataframes
df = dict_to_df(psd_d)
print('Exporting PSD csv')
df.to_csv(outfile, index=False)
sdf = dict_to_df(sig_d)
print('Exporting raw signal csv')
sdf.to_csv(outsig, index=False)

### Get out while you can
print("YOU ARE ALL FREE NOW")
