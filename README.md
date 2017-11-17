_Disclaimer: no actual lasers, just their data._

Python scripts for handling a specific use case of the JPK force-sensing opical trap/Nanotracker 2.

----
# Notices
- The script appears to work for both optical trap _and_ afm data (on linux at least). It detects this based off of the id string in the force-save file: `# bead-id` (optical) vs `# approachID` (afm). As such, I'm using a single script (`fftz.py`) for development going forward. 

# Current State of the Script
- Inputs:
  - Directory/subdirectories containing 1+ optical trap or afm force-save txt files
- Outputs 3 .csv files:
  1. `*_sig.csv` the filtered baseline corrected time series signals
  2. `*_psd.csv` the full psd of the time series above
  3. `*_psd_peaks.csv` the detected peaks of the psd above
- Peak detection thresholds are static and can be set through command line args (but not if script is run as a clickable executable).

----
# Goals
- Implement peak detection using a moving windowed mean/sd threshold

# Roadmap
1. ~~Copy static mean `detect_peaks` implementation from `fftz_afm.py` to `fftz_trap.py`.~~
3. ~~Include code for `detect_peaks.py` in fftz scripts~~
1. Windows test
2. Rolling window threshold for peak detection
4. Windows test
5. Gui parameter selection? _(someday?)_

----
# Non-essential but would be nice
- A freaking stout log
- When run from Windows explorer, if shit crashes, some way of preserving stout (or even just the error message) instead of the terminal window just dissapearing
- Add headers to output files including things like
  - detected force-save type, Fs, scan length, date, etc
  - filter and peak parameters
  - date run, script version
- When printing paths to stout, make them relative to toplevel experiment directory
