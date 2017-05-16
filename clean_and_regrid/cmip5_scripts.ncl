undef("copy_VarAtts_except")
 procedure copy_VarAtts_except(var_from,var_to, Except [*]:string)
 local att_names, i
 begin
     att_names = getvaratts(var_from);
    if(.not.all(ismissing(att_names)))
        do i = 0,dimsizes(att_names)-1
            if (.not.any(att_names(i).eq.Except)) then
                 if (isatt(var_to,att_names(i))) then
                     delete(var_to@$att_names(i)$) ;var_from att may be d      iff size/type
                 end if
                 var_to@$att_names(i)$ = var_from@$att_names(i)$
             end if
         end do
     end if

end     ; of copy_VarAtts_except




;-------------------------------------------------------------------
undef("uniq")
function uniq(x0[*], sort_order)
  ;
  ; This functions extracts all unique elements in a 1D array
  ; of strings or numeric values
  ;
  ; x : 1D array (string or numeric)
  ; sort_order : -1->decreasing, 1->increasing
  ;            It actually seems to be the other way round (2012-11-19)
  ;            - no effect for string input (use 0)
  ;
  ; written 2011-10-26
local x0, x, sort_order, hold_x, Sx, y, i, idx, hold_y, dummy
begin

  ; as input parameters can not be deleted
  x = x0

  ; missing values lead to problems
  if isnumeric(x) then
    idx    = ind(.not.ismissing(x))
    hold_x = x(idx)
    delete([/idx,x/])
    x      = hold_x
    delete(hold_x)
  end if


  Sx = dimsizes(x)
  y = x(0)
  
  do i = 1,Sx(0)-1
    idx = ind(y.eq.x(i))

    if all(ismissing(idx)) then
      hold_y = array_append_record(y,x(i),0)
      delete(y)
      y = hold_y
      delete(hold_y)
    end if

    delete(idx)
  end do

  ; sort the result
  if isnumeric(x) then
    dummy = dim_pqsort(y,sort_order*-2)
  else
    sqsort(y)
  end if

  return(y)
end

;-------------------------------------------------------------------
undef("mdl_rip_xs")
function mdl_rip_xs(x1[*]:string,x2[*]:string)
  ;
  ; Take the cross section of x1 and x2
  ;
  ; I use this function to get those realizations
  ; that are available for the historical and a
  ; certain RCP experiment, i.e.
  ; mdl_rip_0 = mdl_rip_xs(mdl_rip_rcp,mdl_rip_hist)
  ;
local x1,x2,S_x1,idx,xs
begin

  S_x1 = dimsizes(x1)
  do i = 0,S_x1(0)-1

    idx = ind(x2.eq.x1(i))
    if .not.all(ismissing(idx)) then

      if isdefined("xs") then
        hold_xs = array_append_record(xs,x1(i),0)
        delete(xs)
        xs = hold_xs
        delete(hold_xs)
      else
        xs = x1(i)
      end if
    end if
    delete(idx)
  
  end do

  return(xs)

end

;-------------------------------------------------------------------
undef("exclude_long_run")
function exclude_long_run(files[*]:string)
  ;
  ; Get rid of all the files after the 21st century,
  ; i.e. the "Extended Concentration Pathway"-files
  ; Do this by assuming that no files of the 21st
  ; century end after 2100
  ;
  ; input:
  ; files: string array with cmip5 data file names
  ;
  ; Works under the assumption that at least one
  ; file in files ends before or in 2100
  ;
local files, date_range, enddate_nc, endyear, idx, result
begin

  date_range = str_get_field(files,6,"_")
  enddate_nc = str_get_field(date_range,2,"-")
  endyear    = tointeger(str_get_field(enddate_nc,1,".")) / 100
  idx        = ind(endyear.le.2100)
  result = files(idx)

  return(result)

end

;-------------------------------------------------------------------
undef("exclude_mdls")
function exclude_mdls(mdl_rip_in[*]:string,through_outs[*]:string)
  ;
  ; get rid of certain models in mdl_rip array
  ;
  ; mdl_rip : mdl_rip
  ; through_outs : mdl of whiche all realizations will be
  ;           eliminated from mdl_rip
  ;
local mdl_rip_in, mdl_rip, through_outs, St, i_t, idx, hold_mdl_rip
begin

  mdl_rip = mdl_rip_in

  St = dimsizes(through_outs)
  do i_t = 0,St(0)-1
    idx = ind(str_get_field(mdl_rip,1,"_").ne.through_outs(i_t))
    if .not.all(ismissing(idx)) then
      hold_mdl_rip = mdl_rip(idx)
      delete(mdl_rip)
      mdl_rip = hold_mdl_rip
      delete(hold_mdl_rip)
    end if
    delete(idx)
  end do

  return(mdl_rip)
end

;-------------------------------------------------------------------
undef("get_mdl_colors")
function get_mdl_colors(mdl_rip[*]:string,c1[1]:integer,c2[1]:integer)
  ;
  ; x     : 1D mdl_rip array
  ; c1,c2 : Color range, have to be within the
  ;         the number of colors of colortable
  ;         i.e. max [2,255]
  ;
  ; returns a color array, the size of mdl_rip,
  ; with the same color for all realizations
  ; from one model
  ;
local mdl_rip, c1, c2, S_mdl_rip, colors, c_mdl, \
      c_uniq_mdl, S_c_uniq_mdl, mdl_color, count
begin

  S_mdl_rip  = dimsizes(mdl_rip)
  colors = new(S_mdl_rip(0),integer)

  c_mdl = str_get_field(mdl_rip,1,"_")
  c_uniq_mdl = uniq(c_mdl,0)
  S_c_uniq_mdl = dimsizes(c_uniq_mdl)
  mdl_color = tointeger(fspan(c1,c2,S_c_uniq_mdl(0)))

  count = 0
  do j = 0,S_mdl_rip(0)-1
    do i = 0,S_c_uniq_mdl(0)-1

      if isStrSubset(mdl_rip(j),c_uniq_mdl(i)+"_") then
        colors(count) = mdl_color(i)
        count = count+1
      end if

    end do
  end do

  return(colors)
end

;-------------------------------------------------------------------
undef("get_run1")
function get_run1(mdl_rip[*]:string)
  ;
  ; mdl_rip : mdl_rip will be reduced to one run per model.
  ;           The realization with the lowest number is chosen.
  ;
local mdl_rip, mdl, rip, uniq_mdl, S_uniq_mdl, result, idx, \
      rip_tmp, r0, runs, idx2
begin

  mdl = str_get_field(mdl_rip,1,"_")
  rip = str_get_field(mdl_rip,2,"_")

  uniq_mdl = uniq(mdl,0)
  S_uniq_mdl = dimsizes(uniq_mdl)

  result = uniq_mdl+"_"
   
  do i = 0,S_uniq_mdl(0)-1
    
    idx     = ind(mdl.eq.uniq_mdl(i))
    rip_tmp = rip(idx)

    r0   = str_get_field(rip_tmp,1,"p")
    runs = tointeger(str_get_field(r0,1,"r"))
    idx2 = ind(runs.eq.min(runs))

    result(i) = result(i)+rip_tmp(idx2(0))
    
    delete([/idx,rip_tmp,r0,runs,idx2/])
  end do


  return(result)
end


;-------------------------------------------------------------------
undef("files_in_period")
function files_in_period(files[*]:string,y1,y2)
  ;
  ; files : A list of files with data from one model, one variable,
  ;         from one realization. 
  ; y1    : start year of the interval of interest
  ; y2    : end year of the interval of interest
  ;
  ; Result : Reduces the list of files to those which contain data
  ;          in the period of interest, i.e. [y1,y2]
local files, y1, y2, files_dim, pin, time, opt, y1_time, y2_time, \
      hold_result, result

begin
  
  ; get number of files
  files_dim = dimsizes(files)
  
  ; run through files
  do i_f = 0,files_dim(0)-1

    ; convert y1 and y2 to the time units of the file
    ; take care of "360" calendars
    pin = addfile(files(i_f),"r")
    time = pin->time
    opt = 1
    opt@calendar = time@calendar

    y1_time = cd_inv_calendar(y1,1,1,0,0,0,time@units,opt)
    if isStrSubset(opt@calendar,"360") then
      y2_time = cd_inv_calendar(y2,12,30,23,59,0,time@units,opt)
    else
      y2_time = cd_inv_calendar(y2,12,31,23,59,0,time@units,opt)
    end if

    ; test if there some data of this file are within [y1,y2]
    if (max(time).lt.y1_time) .or. (min(time).gt.y2_time) then
      delete(time)
      continue
    end if
    delete(time)

    ; if so, add file to list of files to be read
    if isdefined("result") then
      hold_result = array_append_record(result,files(i_f),0)
      delete(result)
      result = hold_result
      delete(hold_result)
    else
      result = files(i_f)
    end if

  end do ; /i_f

  ; return the eventually shortend list of files
  return(result)
end


;-------------------------------------------------------------------
undef("same_time_coord")
function same_time_coord(fin[*]:string)
  ;
  ; files : A list of files with data from one model, most likely 
  ; with one variable from one realization. 
  ;
  ; Result : True if the time variable in all files has the same
  ;           units, i.e. "days since 1850-1-1".
  ;          False if the above is not the case, e.g. if the time
  ;          units are "days since the start of the file"
local fin, result, S_fin, pin, units0, i_f, units

begin
  
  result = True
  S_fin = dimsizes(fin)

  ; if files contains only one file
  if S_fin(0).eq.1 then
    return(result)
  else

    ; get the time units of the 1st file
    pin    = addfile(fin(0),"r")
    units0 = pin->time@units

    do i_f = 1,S_fin(0)-1

      ; get the time units of the other files
      pin   = addfile(fin(i_f),"r")
      units = pin->time@units

      ; compare, set to False if not all units are equal
      if units.ne.units0 then
        result = False
        break
      end if 

    end do ; i_f

    return(result)
  end if 
end


;-------------------------------------------------------------------
undef("time_file_units")
function time_file_units(fin[*]:string)
  ;
  ; files : A list of files with data from one model, most likely 
  ; with one variable from one realization. 
  ;
  ; Result : time, with units equals to the units of the
  ;          first file of fin
  ;        
  ;        
local fin, result, S_fin, pin, units0, i_f, units

begin
  
  S_fin = dimsizes(fin)

  ; get the time of the 1st file
  pin    = addfile(fin(0),"r")
  time  = pin->time

  ; get the target meta data
  opt = 1
  opt@calendar = time@calendar
  units        = time@units

  do i_f = 1,S_fin(0)-1

    ; read the next file 
    pin   = addfile(fin(i_f),"r")
    time0 = pin->time

    ; convert to 1st file units
    t     = cd_calendar(time0,-5)
    time1 = cd_inv_calendar(t(:,0),t(:,1),t(:,2),t(:,3),t(:,4),t(:,5),units,opt)
    delete([/time0,t/])

    ; add time together
    hold_time = array_append_record(time,time1,0)
    delete([/time,time1/])
    time = hold_time
    delete(hold_time)

  end do ; i_f

  if iscoord(time,"time") then
    delete(time&time)
  end if

  return(time)
end






