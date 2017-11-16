_Warning: no actual lasers, just their data._

Python codes for handling a specific use case of the JPK force-sensing opical trap/Nanotracker 2.

----
# Notices
- Peak detection functionality is dependent on `detect_peaks.py`, meaning that this file _must be in the same directory from which you launch the script or otherwise be in your python path._ 
- The script appears to work for both optical trap _and_ afm data (on linux at least). It detects this based off of the id string in the force-save file: `# bead-id` (optical) vs `# approachID` (afm). As such, I'm using a single script (`fftz.py`) for development going forward. 

# Goals
- Implement peak detection using a moving windowed mean/sd threshold

# Roadmap
1. ~~Copy static mean `detect_peaks` implementation from `fftz_afm.py` to `fftz_trap.py`.~~
2. Rolling window threshold for peak detection
3. ~~Include code for `detect_peaks.py` in fftz scripts~~
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
