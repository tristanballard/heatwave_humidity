Download data using wget scripts created by CMIP website.

Regrid the data using cmip5to1deg.ncl. This must be in the ESMF_regrid folder on Sherlock.
This runs into memory issues frequently, especially if there are multiple runs available.
For CSIRO with 10 runs, I actually did them 3 runs at a time by transferring all the raw data
out of the folder save for 3 runs, temporarily, b/c the code works by reading in all data
available in the raw data folder. Even with 3 runs at once it needed around 250GB memory.

The tristan.py was a handy script to make a cmip5to1deg.ncl file for each of the cmip models and scenarios, 
instead of doing them one by one. It concatenates the cmip5to1deg.historical.p1.txt and cmip5to1deg.historical.p2.txt files then
makes one for each of the cmip models, appropriately named.

After regridding, convert to .rds file format if reading in the data in R b/c it was being
buggy with such big .ncl files. 



