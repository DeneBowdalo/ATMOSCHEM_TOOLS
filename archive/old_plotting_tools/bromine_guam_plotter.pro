;Dene Bowdalo 15/02/2013
WINDOW, /FREE, XSIZE=1500, YSIZE=800

colour_range =[1,2,3,4,5,6,7]
counter = 0
bromine_species = [44,45,46,47,48,49,50,51,52,53]
bromine_title = ['Br2', 'Br', 'BrO', 'HOBr', 'HBr', 'BrNO2', 'BrNO3', 'CHBr3', 'CH2Br2', 'CH3Br'] 
; Choose choors to iterate for the species
for i = 2006,2012 do begin

   add = STRTRIM(i, 2)
   FileName  = 'ctm.bpch.'+add+'0101'
   
   Latitude       = 13.45 
   Longitude       = 144.78 

;Access gridinfo and modelinfo data

   ;modelinfo = CTM_TYPE('GEOS5',resolution=2)
   ;gridinfo  = CTM_GRID(modelinfo)

;Obtain lat and lon of location interested in within model using actual lat/lon
   ;sorted_lons = sort(abs(gridinfo.xmid-144.78))
   ;sorted_lats = sort(abs(gridinfo.ymid-13.45))

   ;model_lon   = sorted_lons(0)
   ;model_lat   = sorted_lats(0)

   ;model_lon_index = where(gridinfo.xmid eq model_lon)
   ;model_lat_index = where(gridinfo.ymid eq model_lat)


;If wish to iterate over whole set bromine species:
    iterate = 'yes'
    if iterate eq 'yes' then begin
    for bromine_species = 44,53 do begin
    if i eq 2006 then begin
    ctm_plot, 'IJ-AVG-$', Filename = Filename, Tracer = bromine_species[counter], Lon = Longitude, Lat = Latitude,  Altrange = [0,20],AVERAGE=3,  color = colour_range[counter], linestyle = 2, title = 'Br for GEOS-v90103-2x5 over averaged Jans in Guam'
    endif else begin 
    ctm_plot, 'IJ-AVG-$', Filename = Filename, Tracer = '45', Lon = Longitude, Lat = Latitude,  Altrange = [0,20],AVERAGE=3,/overplot,  color = colour_range[counter], linestyle = 2
    endelse
    counter2 = counter2+1§
    endfor
    endif
    

if iterate ne 'yes' then begin    
if i eq 2006 then begin
ctm_plot, 'IJ-AVG-$', Filename = Filename, Tracer = bromine_species[counter], Lon = Longitude, Lat = Latitude,  Altrange = [0,20],AVERAGE=3,  color = colour_range[counter], linestyle = 2, title = 'Br for GEOS-v90103-2x5 over averaged Jans in Guam'
endif else begin 
ctm_plot, 'IJ-AVG-$', Filename = Filename, Tracer = '45', Lon = Longitude, Lat = Latitude,  Altrange = [0,20],AVERAGE=3,/overplot,  color = colour_range[counter], linestyle = 2
endelse
counter = counter+1
endfor
endif

LEGEND,         LABEL=['2006','2007','2008', '2009', '2010','2011','2012'], $
                HALIGN=0.95,                            $
                VALIGN=0.95,                           $
                FRAME =1,                              $
                CHARSIZE=1.2,                          $
                LINE = [2,2,2,2,2,2,2],                $
                LCOLOR=[1,2,3,4,5,6,7]
end

