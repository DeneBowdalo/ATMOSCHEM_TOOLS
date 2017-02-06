;Dene Bowdalo 15/02/2013
;WINDOW, /FREE, XSIZE=1500, YSIZE=800

colour_range =[1,2,3,4,5,6,7]
counter = 0
bromine_species = [44,45,46,47,48,49,50,51,52,53]
bromine_title = ['Br2', 'Br', 'BrO', 'HOBr', 'HBr', 'BrNO2', 'BrNO3', 'CHBr3', 'CH2Br2', 'CH3Br'] 
species_h_leg_pos = [0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.05,0.05]
species_v_leg_pos = [0.95,0.05,0.05,0.13,0.95,0.05,0.05,0.95,0.05,0.05]

location = 'Nauru'

; Choose choors to iterate for the species
counter2=0
for bromine_species = 44,53 do begin
set_plot,'ps'
device,filename='CAST_bromine_plots/'+location+'/'+bromine_title[counter2]
device, /color
device, bit=8

counter = 0

for i = 2006,2012 do begin

   add = STRTRIM(i, 2)
   FileName  = 'ctm.bpch.'+add+'0101'
   
   Latitude       = -0.5273
   Longitude       = 166.9367 



if i eq 2006 then begin
ctm_plot, 'IJ-AVG-$', Filename = Filename, Tracer = bromine_species, Lon = Longitude, Lat = Latitude,  Altrange = [0,20],AVERAGE=3,  color = colour_range[counter], linestyle = 2, title =bromine_title[counter2] + ' for GEOS-v90103-2x5 over averaged Jans in ' +location, /NOERASE
endif else begin 
ctm_plot, 'IJ-AVG-$', Filename = Filename, Tracer = bromine_species, Lon = Longitude, Lat = Latitude,  Altrange = [0,20],AVERAGE=3,/overplot, /NOERASE,  color = colour_range[counter], linestyle = 2
endelse
counter = counter + 1
endfor    
LEGEND,         LABEL=['2006','2007','2008', '2009', '2010','2011','2012'], $
                HALIGN=species_h_leg_pos[counter2],                            $
                VALIGN=species_v_leg_pos[counter2],                           $
                FRAME =1,                              $
                CHARSIZE=1.2,                          $
                LINE = [2,2,2,2,2,2,2],                $
                LCOLOR=[1,2,3,4,5,6,7]

device,/close

counter2 = counter2+1
endfor
end

