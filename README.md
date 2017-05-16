# heatwave_humidity

This contains the general scripts for analyzing specific humidity and temperature maximum during heat waves, their trends, and their joint probability. The reanalysis daily data are freely available from NCEP/DOE Reanalysis II. CMIP5 historical (~1950-2005) and CMIP5 rcp85 (2006-2100) data are freely available from CMIP portals. 

To download the CMIP5 data use the wget scripts or download manually.
After downloading CMIP5, regrid all onto a common 1x1 degree grid using the scripts in 'clean_and_regrid'.

'station_data' looks at station observations for temperature and humidity.

Most R scripts are run on a cluster and contain in begginging comments approximate run time and required memory allotment. They are often run in parallel, especially for the CMIP5 data.

Home/data directories will need to be changed accordingly.
