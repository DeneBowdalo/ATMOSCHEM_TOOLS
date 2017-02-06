;read in 2 files
;set_plot,'ps'
;device,filename='/home/db876/model_older/23_11_12_longrun/run/sulphurs_long.ps'
;device, /color
;device, bit=8
;multipanel, 4

;files= file_search('/home/db876/Model_v90103/21_12_12_2x25_long/run/ctm.bpch.2*')


;for i=2006,2006 do begin
 ;   for x=2,2 do begin
;l =string(i)
;k =string(x,format='(i02)')
;date=string(l+k+'01',format='(i08)')
;ctm_get_data, model1, file='/home/db876/Model_v90103/21_12_12_2x25_long/run/ctm.bpch.'+date
 ;   endfor
;endfor

;ctm_get_data, model1, file='/home/db876/Model_v90103/21_12_12_2x25_long/run/ctm.bpch.20060201'
bromines= [50,51,52,53] 


for i = 0,5 do begin

 

FileName  = 'ctm.bpch.20120101'
   Lat       = [  13.45 ]
   Lon       = [ 144.78 ]
   Lev       = 38
   Success   = CTM_Get_DataBlock( Data, 'IJ-AVG-$',                 $
                           XInd=XMid, YInd=YMid, ZInd=ZMid,         $
                           Use_FileInfo=FileInfo,                   $
                           Use_DataInfo=DataInfo,                   $
                           ThisDataInfo=ThisDataInfo,               $
                           Tracer=bromines(i),           FileName=FileName,   $
                           GridInfo=GridInfo,  ModelInfo=ModelInfo, $
                           Lev=Lev,            Lat=Lat,             $
                           Lon=Lon,            WE=WE,               $
                           SN=SN,              UP=UP )


print, data 

endfor
  ; if ( not Success ) then return, 0


;create a list of tracers to plot out
;bromines= ['Br2','Br','BrO','HOBr']


; look from 0 to 3 which is the number of tracers we are interested in
;for i=0,3 do begin

;find right field in model1 for tacer interested in, in the right
;diagnostic category, and at the rifght time.
;k=where(model1.tracername eq bromines(i) and $
   ;        model1.category eq 'IJ-AVG-$') 
           ;model1.tau0 eq 184080.00)     

;l=where(model2.tracername eq sulphurs(i) and $
          ; model2.category eq 'IJ-AVG-$' and $
           ;model2.tau0 eq 219144.00)

;cut out model1 and model2 structure fields to plot.
   ;field1=*model1(k(0)).data
   ;field2=*model2(l(0)).data
   
;calculate difference between 2 fields
;diff=field1-field2

;print out info about mean rtio of bottom 10 levels of model
;print, sulphurs(i), 'Mean Ratio', mean(field1(*,*,0:9))  
   
 
;cut out the lowest level of the difference  

;latitude = where(field1 eq 14.00)
;longitude = where(field1 eq 145.00)

;print, latitude

;print, longitude

;layer=field1(130,54,0) 

;print, layer

; plot that up and put in a colour bar and plot the continents
;tvmap, layer, /cbar, /continents, /countries, divisons=8 , title=sulphurs(i), xtitle='Latitude', ytitle='Longitude'

;device,/close
;set_plot,'x'

;endfor


 
end









   
