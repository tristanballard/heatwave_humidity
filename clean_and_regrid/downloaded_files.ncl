load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "cmip5_scripts.ncl"


; This routine displays graphically which file are downloaded.
; All this is based on the file names - no contents are tested.
; 
;
; Choose experiment and frequency.
; It is possible to have a look at one model only.
; It is also possible to look at specific variables.
;
;
; Potential bug: with rxixpx, order of px may not be correct


;--------------------------------------------------------
begin

  ; define experiment
  exp_folder  = "historical"
  ;exp_folder  = "rcp85"
  ;exp_folder  = "rcp45"
  ;exp_folder  = "piControl"

  ; decide on frequency
  freq_folder = "monthly"
  ;freq_folder = "daily"
  ;freq_folder = "6hr"

  ; if desired, look at a one model only
  one_model = False
  this_one_model = "CCSM4"

  ; if desired, look at specific variable only
  select_var = True
  my_var = (/"snw","pr","prsn","mrro","snm"/); (/"tas","zg","mrso"/)

  ;--------------------------------------------------------
  ; these two options don't really have to be changed

  ; 2006-2100 vs. 2101-2300
  long_run = False

  ; regridded or not
  ;version = "original"
  version = "1deg"

  ;--------------------------------------------------------

  if exp_folder.eq."historical" then
    y1_0 = 1850
    y2_0 = 2005
  end if

  if exp_folder.eq."piControl" then
    y1_0 = 0
    y2_0 = 5000
  end if
  
  if isStrSubset(exp_folder,"rcp") then
    if long_run then
      y1_0 = 2000
      y2_0 = 2300
    else
      y1_0 = 2006
      y2_0 = 2100
    end if
  end if

  span = y2_0-y1_0
  ;--------------------------------------------------------

  path = "$PI_SCRATCH/CMIP5/"+exp_folder+"/"+freq_folder+"/"+version+"/"
  years = ispan(y1_0,y2_0,1)
  Sy = dimsizes(years)
  ny = Sy(0)

  ; get files and strip away the path
  files = systemfunc("ls "+path+"*.nc")
  nd = str_fields_count(files,"/")
  files = str_get_field(files,nd(0),"/")
  Sfiles = dimsizes(files)
  nf = Sfiles(0)

  ; decode the file name
  variable   = str_get_field(files,1,"_")
  realm      = str_get_field(files,2,"_")
  frequency  = str_get_field(files,2,"_")
  model      = str_get_field(files,3,"_")
  experiment = str_get_field(files,4,"_")
  rip        = str_get_field(files,5,"_")
  time       = str_get_field(files,6,"_")
  time1      = tochar(str_get_field(time,1,"-"))
  time2      = tochar(str_get_field(time,2,"-"))
  y1m1       = tostring(time1(:,0:5))
  y2m2       = tostring(time2(:,0:5))

  ; decode rip further
  r = tointeger(str_get_field(str_get_field(rip,1,"i"),1,"r"))
  p = tointeger(str_get_field(rip,2,"p"))

  ; extract the available models and variables
  ; use my uniq_st for this
  all_var = uniq(variable,0)
  all_mdl = uniq(model,0)
  all_rip = uniq(rip,0)

  
  ; if desired, confine to one model
  if one_model then
    
    ; get the available rip's
    f_1m        = systemfunc("ls "+path+"*"+this_one_model+"*.nc")
    rip_1m = uniq(str_get_field(f_1m,5,"_"),0)
    S_rip_1m = dimsizes(rip_1m)

    ; reorder rip in case there is a r10
    tmp = str_get_field(rip_1m,1,"r")
    tmp = str_get_field(tmp,1,"i")
    r10_idx = dim_pqsort(tointeger(tmp),1)
    rip_1m = rip_1m(r10_idx)
    delete([/tmp,r10_idx/])
    
    ; confine to model
    idx = ind(all_mdl.eq.this_one_model)
    hold = all_mdl(idx)
    delete(all_mdl)
    all_mdl = hold
    delete([/idx,hold/])
  end if

  ; if desired, confine to predefined variables
  if select_var then
    delete(all_var)
    all_var = my_var
  end if

  ; 
  S_all_var = dimsizes(all_var)
  S_all_mdl = dimsizes(all_mdl)
  S_all_rip = dimsizes(all_rip)
  n_var     = S_all_var(0)
  n_mdl     = S_all_mdl(0)
  n_rip     = S_all_rip(0)

  ; coverage array
  x = new((/n_var,n_mdl,n_rip,ny/),float)
  x!0 = "variable"
  x!1 = "model"	
  x!2 = "realization"
  x!3 = "year"
  x&year = years


  ; fill the coverage array
  do i = 0,nf-1

  
    idx1 = ind(variable(i).eq.all_var)
    idx2 = ind(model(i).eq.all_mdl)
    idx3 = ind(rip(i).eq.all_rip)
 
    if any(ismissing((/idx1,idx2,idx3/))) then
      continue
    end if

    y1 = tointeger(y1m1(i)) / 100
    y2 = tointeger(y2m2(i)) / 100

    ; skip if necessary
    if y1.gt.y2_0 then
      continue
    end if

    ; and shorten if necessary
    y2 = min((/y2,y2_0/))

    ; assign the y values
    x(idx1,idx2,idx3,{y1:y2}) = (idx2+1)*100 - (r(i)-1)*8 - (p(i)-1)*2

  end do

  ; coverage array 4D -> 2D
  x1d = ndtooned(x(model|:,realization|:,variable|:,year|:))
  x3d = onedtond(x1d,(/n_mdl*n_rip,n_var,ny/))

  x1d = ndtooned(x3d)
  x2d = onedtond(x1d,(/n_mdl*n_rip,n_var*ny/))

  ;---------------------------------------------------------
  ; make plot
  wks_file = exp_folder+"_"+freq_folder+"_"+version+"_"+ \ 
             systemfunc("date "+str_get_dq()+"+%Y-%m-%d"+str_get_dq())
  wks      = gsn_open_wks("newpdf",wks_file)

  gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

  res = True
  res@gsnDraw = False
  res@gsnFrame = False

  ; Tickmarks
  res@tmYROn = False
  res@tmXTOn = False
  res@tmYLMinorOn = False
  res@tmXBMinorOn = False

  res@tmXBMode   = "Explicit"
  res@tmXBValues = ispan(0,n_var-1,1)*ny + floor(ny/2.)
  res@tmXBLabels = all_var

  res@tmYLMode   = "Explicit"
  if one_model then
    x_tmp  = ndtooned(x2d)
    idx    = ind(.not.ismissing(x_tmp))
    values = uniq(x_tmp(idx),1)
    res@tmYLValues = values
    res@tmYLLabels = rip_1m
  else
    res@tmYLValues = ispan(1,n_mdl,1)*100
    res@tmYLLabels = all_mdl
  end if

  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01
  res@tmXBLabelAngleF      = 60

  res@tiMainString = exp_folder +" - "+ freq_folder

  res@xyMonoLineColor  = False
  res@xyDashPattern    = 0
  res@xyLineThicknessF = 2.

  res@trXMinF = 0
  res@trXMaxF = n_var * span + 20

  ; NCL constrains limit the number of lines on one page to
  ; 256. Handle therefore when  n_mdl*n_rip>256
  n_pages = (n_mdl*n_rip / 256) + 1

  ; page index 
  p_idx = new((/n_pages,2/),integer)
  p_idx(:,0) = 0
  p_idx(:,1) = n_mdl - 1

  if n_pages.gt.1 then

    mdl_per_page = n_mdl / n_pages

    do i_p = 0, n_pages-2
      p_idx(i_p,0) = i_p*mdl_per_page
      p_idx(i_p,1) = (i_p+1)*mdl_per_page - 1
    end do

    p_idx(n_pages-1,0) = i_p*mdl_per_page
  end if

  ; Make plots
  do j = 0,n_pages-1

    ; get indices to split x2d
    i1 = p_idx(j,0)*n_rip
    i2 = p_idx(j,1)*n_rip + n_rip - 1
    n_mdl_j = p_idx(j,1) - p_idx(j,0) + 1


    ; get colors
    colors = new((/n_mdl_j,n_rip/),integer)
    do i = 0,n_rip-1
      colors(:,i) = tointeger(floor(fspan(50,250,n_mdl_j)))
    end do
    res@xyLineColors = ndtooned(colors)

    ; in case of one model I cannot use the color above
    ; as it would be white
    if one_model then
      colors_one = ndtooned(colors)
      colors_one = 104
      res@xyLineColors = colors_one
      delete(colors_one)
    end if

    ; make nice range
    res@trYMinF = min(x2d(i1:i2,:))-20
    res@trYMaxF = max(x2d(i1:i2,:))+10
    ; make plot
    plot = gsn_csm_xy(wks,ispan(1,n_var*ny,1)*1.,x2d(i1:i2,:),res)

    ; add vertical lines
    dummy = new(n_var,graphic)
    do i = 0,n_var-1
      dummy(i) = gsn_add_polyline(wks,plot,(/ny*(i+1),ny*(i+1)/),(/0,n_mdl*100/),False)
    end do

    ; plot panels
    gsn_panel(wks,plot,(/1,1/),False)

    delete(res@xyLineColors)
    delete([/dummy,colors/])
  end do


end
