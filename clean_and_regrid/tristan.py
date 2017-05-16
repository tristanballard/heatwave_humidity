import os

List = ["ACCESS1-3","ACCESS1-0","CCSM4","CNRM-CM5","CSIRO-Mk3-6-0","CanESM2","FGOALS-g2","FGOALS-s2","EC-EARTH","GFDL-ESM2M","GFDL-ESM2G","GISS-E2-H","GISS-E2-H-CC","GISS-E2-R-CC","GISS-E2-R","HadCM3","HadGEM2-CC","HadGEM2-ES","IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC4h","MIROC-ESM","MIROC-ESM-CHEM","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P","MRI-CGCM3","NorESM1-M","NorESM1-ME","bcc-csm1-1","bcc-csm1-1-m","inmcm4","BNU-ESM","CMCC-CM","CMCC-CMS","GFDL-CM3","CanCM4","HadGEM2-AO","CESM1-FASTCHEM","CESM1-CAM5","CESM1-BGC","CESM1-WACCM","FIO-ESM"]

submitfile = 'cmip5to1deg.historical.sh'
model_basefile = 'cmip5to1deg.historical.p2.txt'
header = 'cmip5to1deg.historical.p1.txt'


for model in List:
	f = open('vars.txt', 'wb')
	string = '\nmodel_choice = ' + model + '\n'
	f.write(bytes(string, 'UTF-8'))
	f.close()

	newfile = 'cmip5to1deg.' + model + '.historical.ncl'

	newsubmitfile = 'cmip5to1deg.' + model + '.historical.sh'
	os.system('cp ' + submitfile + ' ' newsubmitfile)
	g = open(newsubmitfile, 'a')
	string = 'ncl ' + newfile
	g.write(string)
	g.close()

	os.system('cat ' + header + ' vars.txt ' + model_basefile + ' > ' + newfile)
	os.system('sbatch ' + newsubmitfile)
	os.system('rm vars.txt')
