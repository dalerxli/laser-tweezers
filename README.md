_Disclaimer: no actual lasers, just their data._

Python scripts for handling a specific use case of the JPK force-sensing optical trap/Nanotracker 2.

----
# Current State
- 2 scripts:
  1. `fftz.py`
    - calculate PSD from force-save txt files (AFM or optical data)
  2. `detect_windowed_peaks.py`
    - find peaks in PSD files outputed by `fftz.py`

----
# Roadmap
- Add logging functions to `detect_windowed_peaks.py`

----
# Maybe-Someday Goals
- **Some sort of usage manual**
- When printing paths to stout, make them relative to toplevel experiment directory
- Gui, or even console-based, parameter selection
