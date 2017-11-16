_Warning: no actual lasers, just their data._

Python codes for handling a specific use case of the JPK force-sensing opical trap/Nanotracker 2.

----
# Notices
- Peak detection functionality is dependent on `detect_peaks.py`, meaning that this file _must be in the same directory from which you launch the script or otherwise be in your python path._ 

# Goals
- Implement peak detection using a moving windowed mean/sd threshold

# Roadmap
1. Copy static mean `detect_peaks` implementation from `fftz_afm.py` to `fftz_trap.py`.
2. Rolling window threshold for peak detection
3. Include code for `detect_peaks.py` in fftz scripts
4. Windows test
5. Gui parameter selection? _(someday?)_

----
# Non-essential but would be nice
- A freaking stout log
- When run from Windows explorer, if shit crashes, some way of preserving stout (or even just the error message) instead of the terminal window just dissapearing
- When printing paths to stout, make them relative to toplevel experiment directory
