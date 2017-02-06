;set_plot,'ps'
;device,filename='/home/db876/model_older/23_11_12_longrun/run/sulphurs_long.ps'
;device, /color
;device, bit=8
;multipanel, 4


ctm_get_data, model, file='ctm.bpch.20060101'         
          
        
 
;find right field in model1 for tacer interested in, in the right
;diagnostic category, and at the rifght time.

bromine=where(model.tracername eq 'Br' and $
           model.category eq 'IJ-AVG-$')

;cut out model1 and model2 structure fields to plot.
field1=*model(bromine(0)).data

;print out info about mean rtio of bottom 10 levels of model
;print, sulphurs(i), 'Mean Ratio', mean(field1(*,*,0:9))  


;cut out the lowest level of the difference  

layer=field1(*,*,0) 


LON = 130. + findgen(6)*2.5
LAT = 10. + findgen(6)*2



;plot that up and put in a colour bar and plot the continents
tvmap, layer, LON, LAT, /cbar, /continents, /countries, divisons=8 , title='Br', xtitle='Latitude', ytitle='Longitude', LIMIT=[10,130,30,150]


end
