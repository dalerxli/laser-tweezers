_Disclaimer: no actual lasers, just their data._

Python scripts for handling a specific use case of the JPK force-sensing optical trap/Nanotracker 2.

----
# Current State
- Power Spectral Density (PSD) Scripts:
  1. `fftz.py`
    - Calculate PSD from force-save txt files (AFM or optical data)
  2. `detect_windowed_peaks.py`
    - Find peaks in PSD files outputed by `fftz.py`
- RMS Scripts:
  - 2 exist so to allow an "easy" comparison between filtered/unfiltered
	rms. However, the scripts are identical except for their default
	arguments. Both can be made to behave like the other by passing the
	command line options `--filter` or `--nofilter`.
	```python
	# rms.py
	parser.set_defaults(filter_on=False) 

	# rms_filtered.py 
	parser.set_defaults(filter_on=True) 
	```
  1. `rms_filtered.py`
	- Calculate the rms of highpass filtered (1 Hz low cutoff) signal. This
	  filter is the same one used in the PSD scripts above. This is the
	  recommended one to use.
  2. `rms.py`
    - Calculate the rms of signal without filter application. Could really
	  be called `rms_unfiltered` because the filtered script should be the
	  default use case. Here to allow comparison between
	  filtered/unfiltered rms values.

----
# Planned Additions
- Popup error window for PSD scripts.

----
# Maybe-Someday Goals
- **Some sort of usage manual**
- When printing paths to stout, make them relative to toplevel experiment directory
- Gui, or even console-based, parameter selection
