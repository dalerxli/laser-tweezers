#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python 2.7
""" Calc rms for xyz signals from force-save.txt files, walks through folders
"""


### SCRIPT INFO
__author__ = 'Nick Chahley, https://github.com/pantsthecat/laser-tweezers'
__version__ = '0.1.2'
__day__ = '2018-05-22'
__codename__ = 'Nina Bonita Brown' 

### Import Modules
import csv
import numpy as np
from glob import glob
import os
import sys
import pandas as pd
import datetime
import collections
import logging
import argparse
import scipy.signal

### Arguments 
# I think defining arguments like this is still readable (if a bit cumbersome)
# and allows for a future in which scripts are called from the commandline (in
# a pipeline, or by .exe launchers), even though we are not doing this atm. 
# --Nick
parser = argparse.ArgumentParser()

## Filter toggle and cutoff 
parser.add_argument('-cl', '--cf_low', type=float, default=1.,
    help='Low cutoff freq for bandpass filter in Hz (def 1Hz)')
parser.add_argument('--nofilter', dest='filter_on', action='store_false',
                    help='Highpass filter off')
parser.add_argument('--filter', dest='filter_on', action='store_true',
                    help='Highpass filter on (default)')
parser.set_defaults(filter_on=True) 
parser.add_argument('--noblc', dest='baseline_correct', action='store_false',
                    help='Do not apply baseline correction to signal before rms calculation')

# Access an arg value by the syntax 'args.<argument_name>'
args = parser.parse_args()


### Function Defs
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
def get_cwd():
    ### Check OS switch, this is a bad name b/c of similarity to os.getcwd()
    if sys.platform.lower().startswith('win'):
        cdir = os.getcwd().split('\\')[-1] # Windows
    else:
        cdir = os.getcwd().split('/')[-1] # *.nix
    return cdir
def dirsep():
    """Check OS and set appropriate dir separator 
    """ 
    if sys.platform.lower().startswith('win'):
        dirsep = '\\'
    else:
        dirsep = '/'
    return dirsep
def read_forcesave(f, col=2):
    """ Reads NanoTracker force-save.txt
    Depends on: class CommentedFile
    In:
        f: force-save file to be opened
        col: column of f to be extracted, 2=z1, 5=z2
    Returns:
        sig: 1D time series signal
    """
    # TODO how to tell when we've reached end of file / next seg
    # Probs a more robust way to do this using with open() as :
    csv_file = csv.reader(CommentedFile(open(f, "rb")),
                          delimiter=' ')

    try:
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
        sig = map(float, sig)
        return sig
    except StopIteration:
        print("Problem reading %s" %f)
def get_params(f):
    """ Returns: t, fs, date (ndarray, float, string)
    """
    # TODO num segments -> num loops (only rel for old ramps)
    T = float(get_match_val(f, "settings.segment.1.duration:"))
    N = int(get_match_val(f, "settings.segment.1.num-points:"))
    fs = N/T
    dt = 1/fs
    t = np.arange(dt, T+dt, dt)
    date = get_match_val(f, "bead-id:")
    date = date.split('-')[0] # date is at front of bead-id string
    return t, fs, date

def path_dialog_no_p(whatyouwant):
    """ Prompt user to select a dir (def) or file, return its path
    In
    ---
    whatyouwant : str opts=['folder', 'file']
    """
    import Tkinter
    root = Tkinter.Tk()
    root.withdraw()

    opt = {}
    opt['parent'] = root
    opt['initialdir'] = './'

    if whatyouwant == 'folder':
        from tkFileDialog import askdirectory
        ask_fun = askdirectory
        # dirpath will be to dir that user IS IN when they click confirm
        opt['title'] = 'Please select your experiment directory (be IN this folder)'

    if whatyouwant == 'file':
        from tkFileDialog import askopenfilename
        ask_fun = askopenfilename
        opt['title'] = 'Select psd file to detect peaks from'
        opt['filetypes'] = (('CSV files', '*.csv'), ('All files', '*.*'))

    path = ask_fun(**opt)
    return path

#-------------------------------------------------
### Imported Functions from fftz
## - General logging, param storage/find and headers
# Paramater and Logger functions are all intertwined and interdependent, sry
def path_dialog(whatyouwant):
    """ Prompt user to select a dir (def) or file, return its path

    In
    ---
    whatyouwant : str opts=['folder', 'file']

    Out
    ---
    path : str, Absolute path to file

    """
    import Tkinter
    root = Tkinter.Tk()
    root.withdraw()

    opt = {}
    opt['parent'] = root
    opt['initialdir'] = './'

    if whatyouwant == 'folder':
        from tkFileDialog import askdirectory
        ask_fun = askdirectory
        # dirpath will be to dir that user IS IN when they click confirm
        opt['title'] = 'Please select your experiment directory (be IN this folder)'

    if whatyouwant == 'file':
        from tkFileDialog import askopenfilename
        ask_fun = askopenfilename
        opt['title'] = 'Select psd file to detect peaks from'
        opt['filetypes'] = (('CSV files', '*.csv'), ('All files', '*.*'))

    path = ask_fun(**opt)
    return path
def header_get_params(f, p, kind='rms'):
    """ Where params is a dict containing info/settings 
    
    Issues:
    - date could be different across forcesave files
    """
    
    logger.info("Read params from file %s" %f)

    # TODO exception handling for TypeError
    T = float(get_match_val(f, "settings.segment.1.duration:")) # (s)
    N = int(get_match_val(f, "settings.segment.1.num-points:"))
    fs = N/T
    dt = 1/fs
    ## please don't let this appear on r/shittyprogramming
    p['T'] = T
    p['N'] = N
    p['fs'] = fs

    ## Don't need to calculate these things if rms, but is it even worthwile 
    # to make a diff case for rms to *not* calc them?
    if kind == 'psd':
        p['dt'] = dt
        p['t'] = np.arange(dt, T+dt, dt) # generate time-points for signal (s)
        p['freq'] = freq_calc(len(p['t']), p['fs'])

    id_list = ['bead-id:', 'approachID:']
    date = get_match_val(f, id_list)
    # date is at front of id string ~ "yyyy.mm.dd-hh.mm.ss-*"
    p['date'] = date.split('-')[0] 

    return p
def walk_get_params(path, p):
    """ Read scan params from first force-save encountered and break, assumes 
    all files have the same parameters
    This is a quick fix, and probably bad form but should do alright in a pinch

    "2": uses dict p to store vars
    """
    success = False
    for subdir, dirs, files in os.walk(path):
        if success == False:
            os.chdir(subdir)
            infile = glob('force-save*.txt')
            if len(infile) > 0:
                logger.info('Reading Scan Parameters from dir: %s' %subdir)
                p = header_get_params(infile[0], p, kind='rms')
                success == True
                break
    return p
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
                logger.info('Infering force-save type from dir: %s\n file: %s'\
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

    true_type = [t for t in fs_type.keys() if fs_type[t] == True]
    ## Overbearing mother sanity check
    if len(true_type) == 1:
        logger.info('Found force-save type: %s' %true_type[0])
    else:
        logger.warning('Found multiple force-save types?? : %s' %true_type)

    return fs_type
def premain_setup():
    """ Define params dict p, and setup logger
    """
    # dict to store settings, params and pass 'em to others
    # global p
    p = {}
    p['rootpath'] = path_dialog('folder') # user selects starting folder (QT)

    logger = logger_setup(p) # need rootpath to set logfilename
    logger.info("Version %s (%s) -- \"%s\"" \
                 %(__version__, __day__, __codename__) )
    now = datetime.datetime.now()
    daterun = now.strftime('%Y-%m-%d')
    logger.info('Today is %s' % daterun)
    logger.info('Path to selected experiment folder: %s' %p['rootpath'])

    return p, logger
def logger_setup(p):
    """ Basic setup for crash logger
    """
    logfile = '/'.join((p['rootpath'], os.path.basename(__file__) ))
    logfile = logfile + '.log'
    datefmt = '%H:%M:%S'
    logfmt = '%(asctime)s %(levelname)-8s %(message)s'
    logging.basicConfig(filename=logfile, level=logging.DEBUG,
                        filemode='w', # overwrite log file if exists
                        format=logfmt, datefmt=datefmt)

    ## Console Handler 
    ## Have logger print to stdout as well as log file
    ch = logging.StreamHandler(sys.stdout)
    chfmt = logging.Formatter('%(asctime)s: %(message)s', datefmt)
    ch.setFormatter(chfmt)
    # logger.getLogger().addHandler(ch)
    
    ## !! new testing
    ## boilerplate, "allows per module configuration" 
    logger = logging.getLogger(__name__)
    logger.addHandler(ch)

    return logger
## - Highpass Filter
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

### New Funs for v 0.1x
def dict_to_df_rms(d, cnames):
    df = pd.DataFrame(d.items(), columns = cnames)
    df = df.sort_values(cnames[0]) # sort by scan name
    return df
def scan_transformation_rms(infile, scan_name, output_dict, col=2, filter_hp=False):
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
            # New baseline correction location before filters, if applied
            sig = sig - np.mean(sig)

            ## Run highpass filter, if asked
            if filter_hp == True:
                sig = butter_highpass_filter(sig, highcut=args.cf_low, fs=p['fs'], order=3)
            rms = rms_calc(sig)
            output_dict[colnames[i]] = rms
            logger.info(' '.join(("Finished:", colnames[i])))

        return output_dict
def popup_message(text = '', title='Message'):
    try:
        # Python 3.x imports
        import tkinter as tk
        from tkinter import messagebox
    except ImportError:
        # Fall back to 2.x
        import Tkinter as tk
        import tkMessageBox as messagebox

    root = tk.Tk().withdraw()  # hide the root window

    messagebox.showinfo(title, text)  # show the messagebox
def rms_calc(sig):
    rms = np.sqrt(np.mean(abs(sig)**2))
    return rms
def main_rms_run_basic(p, filter_hp=False):
    rms_d = {}

    # rootpath already defined in premain_setup()
    p = walk_get_params(p['rootpath'], p)

    # report to logger b/c who knows, a crash could be anywhere, anything
    logger.info('Highpass filter : %s' %filter_hp)

    # Assume the first force save file encountered is the same type (optical
    # trap / afm) as all other files involved in this run.
    # Should do this step in walk_get_params above, but w/e, legacy
    p['forcesave_type'] = detect_forcesave_type(p['rootpath'])

    ## Get the rms
    for subdir, dirs, files in os.walk(p['rootpath']):
        os.chdir(subdir)
        scan_name = get_cwd()
        logger.info(' '.join(("Entering", scan_name)))
        rms_d = scan_transformation_rms(glob('force-save*.txt'), scan_name,
                                        rms_d, filter_hp=filter_hp)

    ## Output csv
    # We want to note on outfile (and output csv?) if filter was used
    dv = 'RMS_filtered' if filter_hp == True else 'RMS'

    # TODO probs better if df has cols something like ['scan_name',
    # 'scan_number', 'rms', 'mean_rms']
    df = dict_to_df_rms(rms_d, cnames=['Scan_Name', dv])
    rootpathname = '/'.join((p['rootpath'], p['rootpath'].split('/')[-1]))
    outfile = '_'.join((rootpathname, dv.lower() + '.csv'))
    logger.info('Exporting %s csv' %dv)
    df.to_csv(outfile, index=False)


def main(p):
    """ Main function to log in case of crash. 
    Not an example elegant programming.
    """
    main_rms_run_basic(p, filter_hp=args.filter_on)

#-------------------------------------------------

## log in case main function crashes
if __name__ == "__main__":
    p, logger = premain_setup()
    try:
        main(p)
        logger.info("YOU ARE ALL FREE NOW")

    except Exception as e:
        logger.exception("BAD END. Main crashed. Error: %s", e)
        popup_message("Main crashed. Error: %s" %e, title='BAD END')

