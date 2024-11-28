# Fitting Fermi data

A wrapper for fitting the fermi data with fermipy. Fermipy needs to be installed, and the latest catalog and spacecraft data should be downloaded. 

The `.env` file specifies all relevant paths. Customize as needed. 

In this setup, the skript takes a table (here called `mastertable.txt`) with the relevant entries (names and coordinates fields as `4FGLName`, `RA, `DEC`). 

Run `get_data.py` first to download the most recent data from the Fermi archive for all entries in your table.

The run_fit.sh skript takes the index of the source in the table as argument. You can also specify if the source should be fitted as a power-law or a log-parabola. 
