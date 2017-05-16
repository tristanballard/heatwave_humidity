load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "cmip5_scripts.ncl"
load "ESMF_regridding.ncl"
;-------------------------------------------------------------
;
; -Regrids the cmip5 data onto a 1 degree lat/lon grid
; -Uses ESMF regridding software
; -Three interpolation options: "bilinear", "patch", or "conserve"
; -Read about them here: http://www.ncl.ucar.edu/Applications/ESMF.shtml
;
; This routine maintains the file nomenclature.
;
; The rouinte can regrid all data in a folder or alterntively
; all data from one model or all data from one variable, or
; both.
;
; 
begin

 ; interpolation_style = "bilinear"
  interpolation_style = "patch"
;  interpolation_style = "conserve"

  yy1 = 2006
  yy2 = 2100 

  one_model = False;True ;Set to 'False' to work with many models
  one_var   = True;False;True

  ; define the experiment
;  exp_folder = "historical"
  ;exp_folder = "rcp45"
  ;exp_folder = "rcp60"
  exp_folder = "rcp85"

  ; define the model for single model regridding (even if all one_model is set to 'False' above, one model here must be uncommented)
  model_choice = "ACCESS1-0"
  ;model_choice = "ACCESS1-3"
  ;model_choice = "CCSM4"
  ;model_choice = "CNRM-CM5"
  ;model_choice = "CSIRO-Mk3-6-0"
  ;model_choice = "CanESM2"
  ;model_choice = "FGOALS-g2"
  ;model_choice = "FGOALS-s2"
  ;model_choice = "EC-EARTH"
  ;model_choice = "GFDL-ESM2M"
  ;model_choice = "GFDL-ESM2G"
  ;model_choice = "GISS-E2-H"
  ;model_choice = "GISS-E2-H-CC"
  ;model_choice = "GISS-E2-R-CC"
  ;model_choice = "GISS-E2-R"
  ;model_choice = "HadCM3"
  ;model_choice = "HadGEM2-CC"
  ;model_choice = "HadGEM2-ES"
  ;model_choice = "IPSL-CM5A-LR"
  ;model_choice = "IPSL-CM5A-MR"
  ;model_choice = "IPSL-CM5B-LR"
  ;model_choice = "MIROC4h"
  ;model_choice = "MIROC-ESM"
  ;model_choice = "MIROC-ESM-CHEM"
  ;model_choice = "MIROC5"
  ;model_choice = "MPI-ESM-LR"
  ;model_choice = "MPI-ESM-MR"
  ;model_choice = "MPI-ESM-P"
  ;model_choice = "MRI-CGCM3"
  ;model_choice = "NorESM1-M"
  ;model_choice = "NorESM1-ME"
  ;model_choice = "bcc-csm1-1"
  ;model_choice = "bcc-csm1-1-m"
  ;model_choice = "inmcm4"
  ;model_choice = "BNU-ESM"
  ;model_choice = "CMCC-CM"
  ;model_choice = "CMCC-CMS"
  ;model_choice = "GFDL-CM3"
  ;model_choice = "CanCM4"
  ;model_choice = "HadGEM2-AO"
  ;model_choice = "CESM1-FASTCHEM"
  ;model_choice = "CESM1-CAM5"
  ;model_choice = "CESM1-BGC"
  ;model_choice = "CESM1-WACCM"
  ;model_choice = "FIO-ESM"

  ;define the variable for single variable regridding
  var_choice = "huss"

  ; define the frequency
  freq_folder = "daily"
  ;-------------------------------------------------------

  setfileoption("nc","Format","LargeFile")
  t0 = stringtoint( systemfunc("date +%s"))
  path_in  = "$PI_SCRATCH/CMIP5/"+exp_folder+"/"+freq_folder+"/original/"
  path_out = "$PI_SCRATCH/CMIP5/"+exp_folder+"/"+freq_folder+"/1deg/"
  
  ; get files and strip away the path
  files = systemfunc("ls "+path_in+"*.nc")
  nd = str_fields_count(files,"/")
  files = str_get_field(files,nd(0),"/")
print(files)

  ; decode the file name
  variable   = str_get_field(files,1,"_")
  realm      = str_get_field(files,2,"_")
  frequency  = str_get_field(files,2,"_")
  model      = str_get_field(files,3,"_")
  experiment = str_get_field(files,4,"_")
  rip        = str_get_field(files,5,"_")
  time       = tochar(str_get_field(files,6,"_"))
  y1m1       = tostring(time(:,0:5))
  y2m2       = tostring(time(:,7:12))
  delete([/files,time/])

  ; decode rip further
  r = tointeger(str_get_field(str_get_field(rip,1,"i"),1,"r"))
  p = tointeger(str_get_field(rip,2,"p"))

  ; extract the available models and variables
  ; use my uniq for this
  all_mdl = uniq(model,0)
  all_var = uniq(variable,0)
  print(all_mdl)
  all_rip = uniq(rip,0)
  
  ; get the number of models, variables and runs
  S_all_mdl = dimsizes(all_mdl)
  S_all_var = dimsizes(all_var)
  S_all_rip = dimsizes(all_rip)
  n_mdl     = S_all_mdl(0)
  n_var     = S_all_var(0)
  n_rip     = S_all_rip(0)

  ;--------------------------------------------------------
  ; get the index for this model
  this_mdl = ind(all_mdl.eq.model_choice)
  
  ;--------------------------------------------------------
  ; get the index for this variable
  this_var = ind(all_var.eq.var_choice)
print(this_var)
  ;--------------------------------------------------------
  ; define latitude of new grid
  nlat          = 180
  lat           = fspan(-89.5,89.5,nlat)
  lat@long_name = "Latitude"
  lat@units     = "degrees_north"

  ; define longitude of new grid
  nlon          = 360
  lon           = fspan(0,359,nlon)
  lon@long_name = "Longitude"
  lon@units     = "degrees_east"

  ; constain loop to one model/variable if desired
  if one_model then
    i_m1 = this_mdl
    i_m2 = this_mdl
    print(""+model_choice)
  else
    i_m1 = 0
    i_m2 = n_mdl-1
  end if

  if one_var then
    i_v1 = this_var
    i_v2 = this_var
    print(""+var_choice)
  else
    i_v1 = 0
    i_v2 = n_var-1
  end if

  do i_m = i_m1,i_m2        ; GCM loop
    do i_v = i_v1,i_v2      ; Variable Loop 
      do i_rip = 0,n_rip-1  ; Realization


        ; get files, if they exist for this i_m/ i_v/ i_rip combination
	    files0 := systemfunc("ls "+path_in+all_var(i_v)+"_*"+all_mdl(i_m)+"_*"+all_rip(i_rip)+"_*.nc")
    

	    if all(ismissing(files0)) then
	        delete(files0)
	    continue
	    end if

        files := files0;files_in_period(files0,yy1,yy2)
        delete(files0)
        print(files)
        
        ; Addfiles
    	Sfiles = dimsizes(files)

        if Sfiles.gt.1 then
          pin := addfiles(files,"r")
          ListSetType(pin,"cat")
          time0 := pin[:]->time
          temp := pin[:]->$all_var(i_v)$ 
        else
          pin   := addfile(files,"r")
          time0 := pin->time
          temp  := pin->$all_var(i_v)$
        end if

            ;---UPDATED BY JUSTIN. THIS SECTIONS IS THE NEW REGRIDDING METHOD FOR NCL >6.2
            Opt                =  True
            Opt@SrcFileName    := "./regrid_files/"+all_mdl(i_m)+"_src_SCRIP.nc"
            Opt@DstFileName    := "./regrid_files/"+all_mdl(i_m)+"_dst_SCRIP.nc"
            Opt@WgtFileName    := "./regrid_files/"+all_mdl(i_m)+"_wgt.nc"
            Opt@ForceOverwrite = True
            Opt@SrcGridMask    := where(.not.ismissing(temp),1,0)
            Opt@InterpMethod   = interpolation_style
            Opt@DstGridType    = "1x1"
            Opt@DstLLCorner    = (/ -89.5d, 0.00d/)
            Opt@DstURCorner    = (/  89.5d, 359d/)
            
            ;-_Run the regrid routine
            x1 = ESMF_regrid(temp,Opt)

            ; change LAT->lat and LON->lon
            x1!0 = "time"
            x1!1 = "lat"
            x1!2 = "lon"
            time = cd_calendar(time0,1)
            x1&time = time
            printVarSummary(x1)


; create new file name by changeing the date
	        nd := str_fields_count(files(0),"/")
	        fout := str_get_field(files(0),nd(0),"/")

            date_old := str_get_field(fout,6,"_")
            d1 := time(minind(time))
            d2 := time(maxind(time))
            date_new := d1+"-"+d2+".nc"
            fout := str_sub_str(fout,date_old,date_new)

	        fout := path_out + fout

            ; open file
	        system("if [ -f "+fout+" ] ; then rm "+fout+" ; fi")
	        pout = addfile(fout,"c")

	        ; set define mode
	        setfileoption(pout,"DefineMode",True)

            ; predefine coordinate variables
            dimNames = (/"time", "lat", "lon"/)
            dimSizes = (/ -1   , nlat , nlon /)
            dimUnlim = (/ True , False, False/)
            filedimdef(pout,dimNames,dimSizes,dimUnlim)

            ; predefine dimensionality of variables
            filevardef(pout, "time" ,typeof(time),"time")
            filevardef(pout, "lat"  ,typeof(lat),"lat")
            filevardef(pout, "lon"  ,typeof(lon),"lon")
            print(getvardims(x1))
            filevardef(pout, all_var(i_v), typeof(x1), getvardims(x1))

            ;copy attributes
            filevarattdef(pout,"time" ,time)
            filevarattdef(pout,"lat"  ,lat)
            filevarattdef(pout,"lon"  ,lon)
            filevarattdef(pout,all_var(i_v),x1)

            ; set write mode
            setfileoption(pout,"DefineMode",False)

            ; write output
            pout->time           = (/time/)
            pout->lat            = (/lat/)
            pout->lon            = (/lon/)
            pout->$all_var(i_v)$ = (/x1/)

            delete([/x1,time0,time/])



          end do ; i_rip  Realizations
    ;      delete([/years,ny/])

        end do ; i_v      Variable
        

    print(all_mdl(i_m))   

      end do ; i_m        GCM

      t1 = stringtoint( systemfunc("date +%s"))
      print("elapsed time: "+((t1-t0)/60.)+" minutes")

      print("done")
    end
