#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python 2.7
from __future__ import division, print_function

# NOTE FOR DEV USE
""" PSD calculator and peak detection for JPK Nanotracker 2 Data 

"""
""" note on future (python 3) division

b/c this was in detect_peaks.py, and we should use non-ambiguous div ops

'//' "floor division" operator (python 2 def for ints)
'/' "true division" operator (python 2 def for floats)

>>> from __future__ import division
>>> 5 // 2 
    2
>>> 5 / 2
    2.5
>>> int(5 / 2)
    2
>>> float(5 // 2)
    2.0
"""


### SCRIPT INFO
__author__ = 'Nick Chahley, https://github.com/pantsthecat/laser-tweezers'
__version__ = '0.1.2'
__day__ = '2017-12-07'
__codename__ = 'Coconut Crab' 
print("Version %s (%s) -- \"%s\"" %(__version__, __day__, __codename__) )

### Imports and defs and arguments {{{
# -------------------------------------------------
import argparse
import csv
import numpy as np
from glob import glob
import os
import sys
import pandas as pd
import scipy.signal

### SCRIPT ARGUMENTS (Lowpass filtering below 1Hz) salim
# These are a bunch of "useful" command line flags/args the utility of which
# will probably never see the light of day b/c (a) everyone wants the run the
# script just by double-clicking it and (b) I have no idea how to program a gui 
# --nikoli
parser = argparse.ArgumentParser()
parser.add_argument('-cl', '--cf_low', type=float, default=1.,
    help='Low cutoff freq for bandpass filter in Hz (def 1Hz)')
parser.add_argument('--nofilter', dest='filter_on', action='store_false')
parser.add_argument('--filter', dest='filter_on', action='store_true',
    help='def')
parser.set_defaults(filter_on=True) #def: exported sig will also be filtered

## Sensor specificity
parser.add_argument('-x', '--xfft', dest='run_x_fft', action='store_true',
    help='Run FFT on X signal')
parser.set_defaults(run_x_fft=False)
parser.add_argument('-y', '--yfft', dest='run_y_fft', action='store_true',
    help='Run FFT on Y signal')
parser.set_defaults(run_y_fft=False)
parser.add_argument('-z', '--zfft', dest='run_z_fft', action='store_true',
    help='Run FFT on Z signal (def)')
parser.add_argument('-nz', '--nozfft', dest='run_z_fft', action='store_false',
    help='DO NOT run FFT on Z signal')
parser.set_defaults(run_z_fft=True)

## Peak detection 
# Can turn it on or off. Please do this. Please don't make two versions, one
# w/ detect_peaks true and one with it false. If you're going to do this at
# least make two executables that call the one script but pass different args.
# These flags just want to be loved -- or even just used...
parser.add_argument('--detect_peaks', dest='detect_peaks', action='store_true',
    help='Output separate file of peaks above <thershold_tbd>')
parser.set_defaults(detect_peaks=True)
# These two args get used as defaults for the function get_peaks
# hence defining these before the functions
parser.add_argument('-tsd', '--thresh_sd', type=float, default=3,
    help='Multiple of signal sd to set min peak amplitude at (def 3)')
parser.add_argument('-mpd', '--peak_space', type=int, default=10,
    help='Value for detect_peaks minimum peak distance (def 10)')

# for 'peak_detection_sep.py'
parser.add_argument('-pf', '--psdfile', type=str, 
    help='Path/to/file. Read psd and detect peaks from this csv file, do not process any force-save files') 

# Access an arg value by the syntax 'args.<argument_name>'
args = parser.parse_args()


## "Original" functions/class of the script
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
# def dirsep():
    # """ This function is no longer used, right?
    # """
    # if sys.platform.lower().startswith('win'):
        # dirsep = '\\'
    # else:
        # dirsep = '/'
    # return dirsep
def get_cwd():
    ### Check OS switch 
    if sys.platform.lower().startswith('win'):
        cdir = os.getcwd().split('\\')[-1] # Windows
    else:
        cdir = os.getcwd().split('/')[-1] # *.nix
    return cdir
def freq_calc(N, fs):    
    """ The line above was commented out
    """
    dt = float(1/fs)
    if N % 2 == 0: 
        ind = np.arange(1,N//2) # future division
    else:
        ind = np.arange(1,N//2+1) # future division  
    # Get corresponding fft freqs 
    freq = np.fft.fftfreq(N, dt)
    # drop dc offset and neg freqs
    freq = freq[ind]
    return freq
def read_forcesave(f, col=2):
    """ Reads AFM force-save.txt << nope, back to old. Now use forcesave-type
    detection functions to call either afm_run() or single_channel_run()

    Depends on: class CommentedFile
    In:
        f: force-save file to be opened
    Returns:
        sig: 1D time series signal
    """
    # Probs a more robust way to do this using with open() as :
    csv_file = csv.reader(CommentedFile(open(f, "rb")),
                          delimiter=' ')
    
    # This loop is to skip the Position Ramp segment
    # Runs until it encounters a blank line
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
    N = len(X)
    
    # Baseline corrected signal
    X = X - np.mean(X)

    # Highpass filter -- rm signals below cutoff freq
    if args.filter_on == True:
        
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
        ind = np.arange(1,N//2) # future division
    else:
        ind = np.arange(1,N//2+1) # future division

    # Get that PSD, son!
    psd = abs(fx[ind])**2    

    # X_output: the filtered (if args.filter_on), baseline corrected time-series
    return psd, X_output
def scan_transformation(infile, scan_name, col=2):
    """ Run the above functions and store psd for each scan in a dict (psd_d)
        Also store extracted raw signal in a dict

        channel/sensor codes: 0,1,2 = x,y,z
    """
    if len(infile) > 0: # if there is at least 1 force-save*.txt
        colnames = []
        infile.sort()
        for i in range(len(infile)): # for each force-save*.txt
            colnames.append(''.join((scan_name, '-', str(i+1))))
            sig = read_forcesave(infile[i], col=col)
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

    # This ifel is so that we can arbitrarily search for both 'bead-id' 
    # (optical trap) and 'approachID' (AFM) -- or more options if need be
    if type(match) == list:
        for line in search.readlines():
            if any(x in line for x in match):
                return line.split()[-1]
                break
    else: # Assume that match is a string
        for line in search.readlines():
            if match in line:
                return line.split()[-1]
                break
def get_params(f):
    """ Returns: t, fs, date (ndarray, float, string)
        Maybe use date for start of file name?
    """
    print("Read params from file %s" %f)

    # TODO exception handling for TypeError
    T = float(get_match_val(f, "settings.segment.1.duration:")) # (s)
    N = int(get_match_val(f, "settings.segment.1.num-points:"))
    fs = N/T
    dt = 1/fs
    t = np.arange(dt, T+dt, dt) # generate time-points for signal (s)

    # NOTE test
    id_list = ['bead-id:', 'approachID:']
    date = get_match_val(f, id_list)
    # # Check b/c AFM files use "approachID:" instead of "bead-id"
    # if date is None: # "bead-id" not found in file
        # date = get_match_val(f, "approachID:")   
    # date is at front of id string ~ "yyyy.mm.dd-hh.mm.ss-*"
    date = date.split('-')[0] 

    return t, fs, date
def butter_highpass(highcut, fs, order=3):
    """ Design a digital highpass Butterworth filter, return filter coefficents
    highcut, fs in Hz
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


## Generalised run functions, script SHOULD be able to accept optical or
# afm force-save files (unverified)
def detect_forcesave_type(path):  
    """ Returns a dict of optical/afm and a bool indicating if the file read
        is of that type. Walks through directory tree, terminates at first 
        successful identification.
    """

    # Use this dict to store bools for afm vs optical trap ft id
    # global fs_type
    fs_type = {'optical' : False,
               'afm' : False}
    fs_type_matches = {'optical' : '# bead-id',
                       'afm' : '# approachID'}

    success = False
    for subdir, dirs, files in os.walk(path):
        if success == False:
            os.chdir(subdir)
            infile = glob('force-save*.txt')
            if len(infile) > 0:
                print('Infering force-save type from dir: %s\n file: %s'\
                      %(subdir,infile[0]))

                search = open(infile[0], 'r')

                for line in search.readlines():
                    for key in fs_type_matches.keys():
                        if line.split(':')[0] == fs_type_matches[key]:
                            fs_type[key] = True
                            success = True
                            break
                    if success == True:
                        break
    return fs_type
def afm_run(sensor='AFM', col=1):
    """ Run over all subdirs. This is a copy of single_channel_run() slightly
        changed in the event of an AFM force-save detection
        Col we're interested in is col 1 (vDisplacement)
    """
    print('\nBeginning an AFM run')
    ### Init data dictionaries
    global psd_d
    psd_d = {}
    psd_d['Hz'] = freq
    global sig_d
    sig_d = {}
    sig_d['Time'] = t
        
    ### Attempted os looping
    for subdir, dirs, files in os.walk(rootpath):
        os.chdir(subdir)
        scan_name = get_cwd()
        print(' '.join(("Entering", scan_name)))
        scan_transformation(glob('force-save*.txt'), scan_name, col=col)

    ### Make and export dataframes
    rootpathname = '/'.join((rootpath, rootpath.split('/')[-1]))
    df = dict_to_df(psd_d)
    outfile = '_'.join((rootpathname, sensor, "psd.csv"))
    print('Exporting PSD csv for %s' %sensor)
    df.to_csv(outfile, index=False)

    sdf = dict_to_df(sig_d)
    outfile = '_'.join((rootpathname, sensor, "sig.csv"))
    print('Exporting raw signal csv for %s' %sensor)
    sdf.to_csv(outfile, index=False)

    # return the psd df for processing by detect_peaks
    return df
def single_channel_run(sensor, col, channel=1):
    """ Run over all subdirs, extracting/transforming signal from <channel>
        Yes, this is inefficient. Fight me. I'm sorry.

        sensor column codes: 0,1,2 = x,y,z
        channel is always 1 right now, but maybe someone'll want 2 also in the
        future 
    """
    print('\nBeginning single run for channel %s' %sensor + str(channel))
    ### Init data dictionaries
    global psd_d
    psd_d = {}
    psd_d['Hz'] = freq
    global sig_d
    sig_d = {}
    sig_d['Time'] = t
        
    ### Attempted os looping
    for subdir, dirs, files in os.walk(rootpath):
        os.chdir(subdir)
        scan_name = get_cwd()
        print(' '.join(("Entering", scan_name)))
        scan_transformation(glob('force-save*.txt'), scan_name, col=col)

    ### Make and export dataframes
    rootpathname = '/'.join((rootpath, rootpath.split('/')[-1]))
    df = dict_to_df(psd_d)
    outfile = '_'.join((rootpathname, sensor, "psd.csv"))
    print('Exporting PSD csv for %s' %sensor)
    df.to_csv(outfile, index=False)

    sdf = dict_to_df(sig_d)
    outfile = '_'.join((rootpathname, sensor, "sig.csv"))
    print('Exporting raw signal csv for %s' %sensor)
    sdf.to_csv(outfile, index=False)

    # return the psd df for processing by detect_peaks
    return df

## Detect Peaks suite
def strip_ext(filename, ret_ext=False):
    """ Strip extension from a filename/path, handle any '.' other than '.ext'
    that might appear in path.
    I'm either lazy, or not smart enough to implement this right at the start.
    """
    filename_strip = '.'.join((filename.split('.')[:-1]))
    ext = filename.split('.')[-1]
    if ret_ext == True:
        return ext
    else:
        return filename_strip
def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):

    """Detect peaks in data based on their amplitude and other features.

    __author__ = "Marcos Duarte, https://github.com/demotu/BMC"
    __version__ = "1.0.4"
    __license__ = "MIT"

    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
        ^ We are not using show/_plot. It's commented out --nikoli
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`
    
    The function can handle NaN's 

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Examples
    --------
    >>> from detect_peaks import detect_peaks
    >>> x = np.random.randn(100)
    >>> x[60:81] = np.nan
    >>> # detect all peaks and plot data
    >>> ind = detect_peaks(x, show=True)
    >>> print(ind)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # set minimum peak height = 0 and minimum peak distance = 20
    >>> detect_peaks(x, mph=0, mpd=20, show=True)

    >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
    >>> # set minimum peak distance = 2
    >>> detect_peaks(x, mpd=2, show=True)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # detection of valleys instead of peaks
    >>> detect_peaks(x, mph=0, mpd=20, valley=True, show=True)

    >>> x = [0, 1, 1, 0, 1, 1, 0]
    >>> # detect both edges
    >>> detect_peaks(x, edge='both', show=True)

    >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
    >>> # set threshold = 2
    >>> detect_peaks(x, threshold = 2, show=True)
    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])

    # Not using show/plot in current fftz.py implementation --nikoli, 16-11-2017
    # if show:
        # if indnan.size:
            # x[indnan] = np.nan
        # if valley:
            # x = -x
        # _plot(x, mph, mpd, threshold, edge, valley, ax, ind)
    return ind
def get_peaks(psd_df, thresh_sd=args.thresh_sd, space=args.peak_space):
    """ Return a df containing only peaks. Minimum peak height is equal to
    'mean + thresh_sd * sd' Adjacent columns of the same freq that do not
    contain peaks will be 0. If no peaks are found at a given freq in any
    columns that row will be omitted from the df.

    Notes:
    - 2 options for filling non-peak cells in the df: zeros or NaN
    - NaN: apparently might allow for easier elimination of "empty" cells, I can
      also see it being easier to count, average, etc the peaks
    - Zeros: going w/ this for now b/c I don't know how much Excel and/or MatLab
      will like importing and encountering 'NaN' while plotting -- which is 
      probably what the primary downstream use of the data from this function
      be. We can always have another function to change 0 >> NaN afterwards --
      god knows we aren't programming for peak effeciency here.
      >>> # Eliminate rows that contain only zeros
      >>> peaks_df = peaks_df.loc[(peaks_df!=0).any(axis=1)]
    """
    # they will probably want detect_peaks to be physically in this file...
    # from detect_peaks import detect_peaks

    # Create a df of zeros in which to store the peaks 
    psd_no_hz = psd_df.iloc[:, range(1, len(psd_df.columns))]
    peaks_df = pd.DataFrame(np.zeros_like(psd_no_hz),
                           index=psd_no_hz.index, columns=psd_no_hz.columns)
    # this df was tmp
    psd_no_hz = ""
    
    # iterate detect_peaks over each psd column in the df and add the peaks to
    # df of zeros at their index
    for i in range(1, len(psd_df.columns)):
        y = psd_df.iloc[:, i]
        if thresh_sd == None:
            # detect all peaks regardless of height and worry about a dynamic
            # thresh comparison later. Also don't worry about using time 
            # effectively, b/c this is slow AF
            thresh=None
            print('Peak detection by space only. This will take a while')
        else:
            # static peak height threshold
            thresh = np.mean(y) + thresh_sd * np.std(y)
        ind = detect_peaks(y, mph=thresh, mpd=space)
        peaks_df.iloc[ind, i-1] = psd_df.iloc[ind, i]
        
    # take only rows in which any column contains a non-zero
    peaks_df = peaks_df.loc[(peaks_df!=0).any(axis=1)]

    # add the Hz col back in
    # peaks_df.index == row names, which are same indexes from psd_df
    tmp = psd_df.iloc[:,0].loc[peaks_df.index]
    # note: pd.DataFrame.insert is automatically an inplace operation
    peaks_df.insert(0, 'Hz', tmp)
    
    return peaks_df
def export_peaks(peaks_df, rootpath, sensor):
    """ Export the peaks df to a csv. Great docstring.

    TODO: if this runs okay, incorporate thresh_sd into the outfile name
    """
    rootpathname = '/'.join((rootpath, rootpath.split('/')[-1]))
    outfile = '_'.join((rootpathname, sensor, "psd_peaks.csv"))
    print('Exporting PSD peaks csv for %s' %sensor)
    print('Exporting to: %s' %outfile)
    peaks_df.to_csv(outfile, index=False)
def export_peaks_from_psdfile(peaks_df, infile_path):
    """ Use if we get psd from a file instead of internally from a fft run. 
    Export the peaks df to a csv. Great docstring.

    TODO: if this runs okay, incorporate thresh_sd into the outfile name
        ^ yeah right
    """
    rootpathname = strip_ext(infile_path)
    outfile = '_'.join((rootpathname, "peaks.csv"))
    print('Exporting PSD peaks csv')
    print('Exporting to: %s' %outfile)
    peaks_df.to_csv(outfile, index=False)


## The one where Ross tries to clean up and compartmentalize different uses 
## for the script but introduces a bunch of global variables instead
def main_fft_run(filter_on = True):
    """In a function so we can conveniently choose to not run it.
    But the detect peaks logic/funs are in it as well so... fuck
    """

    # Report variable settings for fft/psd run
    if filter_on == True:
        print("Butterworth order 3 highpass filter is ON")
        print("Low cutoff frequency is %d Hz" % args.cf_low)
    else:
        print("Butterworth highpass filter is OFF")

    # Backslash is the worst path separation character
    global rootpath
    rootpath = folder_dialog() # user selects starting folder (QT)

    # Open first force-save*.txt file we can find and read/calculate scan 
    # paramaters from the header of that file. *assumption that params are consistant
    # across all scans*
    #    t: time points for time domain signal data (x-axis)
    #    dt: time resolution (1/fs)
    #    freq: the frequency vales (0:Fnyq Hz) for the frequency domain signal
    # if I wasn't a hack I'd do this w/ a dict, not globals
    global t, dt, freq
    t, dt, freq = walk_get_params(rootpath)

    # Building on our house of assumptions: the first force save file encountered
    # is the same type (optical trap / afm) as all other files involved in this run.
    forcesave_type = detect_forcesave_type(rootpath)

    ### Optical Trap Run(s)
    if forcesave_type['optical'] == True:
        print('Detected Forcesave Type: Optical Trap')
        # Very sophisticated logic for deciding which sensors to run
        # Detect peaks logic included w/n the afm/trap logic as a (pointless
        # and ineffective?) attempt at future proofing someone wanting to run
        # multiple trap channels (eg x AND z) in a single execution.
        # Also so we can lamely pass sensor id to export_peaks
        if args.run_x_fft == True:
            sensor = 'x'
            psd_df = single_channel_run(sensor, 0)

        if args.run_y_fft == True:
            sensor = 'y'
            psd_df = single_channel_run(sensor, 1)

        if args.run_z_fft == True:
            sensor = 'z'
            psd_df = single_channel_run(sensor, 2)

    ### AFM Run
    if forcesave_type['afm'] == True:
        print('Detected Forcesave Type: AFM')
        psd_df = afm_run(col=1)
def peaks_from_psdfile(infile_path, space=args.peak_space):
    """ Here we try to make a dynamic windowed threshold for peak height """
    psd_df = pd.read_csv(infile_path, sep=',')
    # get a df of peaks regardless of any height thresholds
    peaks_df = get_peaks(psd_df, thresh_sd=None, space=space)

    # window size should approximate space (+/- 1)
    k = space // 2

    export_peaks_from_psdfile(peaks_df, infile_path)

### Main logic
def main_fftpeaks_logic():
    """ Decide whether to run psd calculation (fft) or peak detection.
    """
    if args.psdfile:
        ## Commenting out atm. Use detect_windowed_peaks.py for 
        ## peak detection
        # peaks_from_psdfile(args.psdfile)
        print('NA. Please use detect_windowed_peaks.py for peak finding')
    else:
        main_fft_run(filter_on = args.filter_on)
#------------------------------------------------- }}}

### NOT FUNCTION DEFINITIONS
# run the regular script
main_fftpeaks_logic()

### Get out while you can
print("YOU ARE ALL FREE NOW")
