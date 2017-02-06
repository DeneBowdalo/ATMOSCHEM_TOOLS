pro readin_gaw_sites, names, locations
openr,1,'sites.dat'
line=''
readf,1,line
names=strarr(100)
locations=fltarr(100,2)
counter=0
while not eof(1) do begin
   readf,1,line
   lister= strsplit(line, string(9b), /extract)
   print, lister(0), float(lister(2)), float(lister(3))
   names(counter)=lister(0)
   locations(counter,*)=float([lister(2:3)])
   counter=counter+1
endwhile
close, 1
names=names(0:counter-1)
locations=locations(0:counter-1,*)
end



;+
; Writes a set of Planefile.dat.YYYYMMDD files for the whole year
; to be read in by GEOS-Chem.v8-01-01.
; 
; Change the year (y), and points to those you desire, i.e. lat, lon,
; pressure, and add a marker in the 'type'.
; 
; Current set up is to output every hour. Check that you have the
; required tracers in 'list'.
;-

itop_reader,'merged-data_faam_201107_r1_b616-b632_v01_060s-all.na',$
  names, x,data, multi,error,notes1, notes2,add

k_lat=where(names eq 'XR5 GPS Latitude')
k_lon=where(names eq 'XR5 GPS Longitude')
k_pre=where(names eq 'Static pressure (hPa)')

slist=[$
      'CH2OOB',$
      'CH3CHOOA',$
      'CH2OO',$
      'CH3CHOO',$
      'O3',$
        'NO2',$
        'NO',$
        'NO3',$
        'N2O5',$
        'HNO4',$
        'HNO3',$
        'HNO2',$
        'PAN',$
        'PPN',$
        'PMN',$
        'R4N2',$
        'H2O2',$
        'MP',$
        'CH2O',$
        'HO2',$
        'OH',$
        'RO2',$
        'MO2',$
        'ETO2',$
        'CO',$
        'C2H6',$
        'C3H8',$
        'PRPE',$
        'ALK4',$
        'ACET',$
        'ALD2',$
        'MEK',$
        'RCHO',$
        'MVK',$
        'SO2',$
        'DMS',$
        'MSA',$
        'SO4',$
        'TRA_1',$
        'TRA_2',$
        'TRA_3',$
        'TRA_4',$
        'TRA_5',$
        'TRA_6',$
        'TRA_7',$
        'TRA_8',$
        'TRA_9',$
         'TRA_10',$
         'TRA_11',$
         'TRA_12',$
         'TRA_13',$
         'TRA_14',$
         'TRA_15',$
         'TRA_16',$
         'TRA_17',$
         'TRA_18',$
         'TRA_19',$
         'TRA_20',$
         'TRA_21',$
         'TRA_22',$
         'TRA_23',$
         'TRA_24',$
         'TRA_25',$
         'TRA_26',$
         'TRA_27',$
         'TRA_28',$
         'TRA_29',$
         'TRA_30',$              
         'TRA_31',$
         'TRA_32',$
         'TRA_33',$
         'TRA_34',$
         'TRA_35',$
         'TRA_36',$
         'TRA_37',$
         'TRA_38',$
         'TRA_39',$
         'TRA_40',$     
         'TRA_41',$
         'TRA_42',$
         'TRA_43',$
         'REA_323',$               ;JNO2
         'REA_324',$            ;JH2O2
         'REA_325',$            ;MP
         'REA_326',$            ;CH2O (RADICAL)
         'REA_327',$            ;CH2O (NON-RADICAL)
         'REA_328',$            ;HNO3
         'REA_329',$            ;HNO2
         'GMAO_TEMP',$          ;
         'GMAO_ABSH',$          ;
         'GMAO_SURF',$          ;
         'GMAO_PSFC',$          ;
         'GMAO_UWND',$          ;
         'GMAO_VWND',$          ;
         'REA_N2O5']
  nvar=n_elements(slist)
; Change the year here:

readin_gaw_sites, sites, locations
lats=locations(*,0)
lons=locations(*,1)
pres=lats*0.+1030.

for yr=2010,2011 do begin
    for m=1,12,1 do begin
        mth=strtrim(m,2)
        if (m lt 10) then begin
            mth='0'+mth
        endif
        for d=01,31,1 do begin
            day=strtrim(d,2)
            if (d lt 10) then begin
                day='0'+day
            endif
            y=strtrim(string(yr,format='(i4)'),2)
            openw,1,'Planeflight.dat.'+y+mth+day
            
            printf,1,'Planeflight.dat -- input file for ND40 diagnostic GEOS_5'
            printf,1,'Helen Macintyre'
            printf,1,'15 Sept 2008'
            printf,1,'-----------------------------------------------'
            printf,1, nvar,'    ! Number of variables to be output',$
              format='(i3,a)'
            printf,1,'-----------------------------------------------'
            
            
            
            for i=0,n_elements(slist)-1 do begin
                printf,1,slist(i)
            endfor
            
            printf,1,'-------------------------------------------------'
            printf,1,'Now give the times and locations of the flight'
            printf,1,'-------------------------------------------------'
            printf,1,'Point  Type DD-MM-YYYY HH:MM     LAT     LON   PRESS'
            
            counter=0
            format='(i5,x,a5,x,i2.2,"-",i2.2,"-",i4,x,i2.2,":",i2.2,x,f7.2,x,f7.2,x,f7.2 )'
     
            
                 
            mn=30
            
            for h=0,23 do begin
                date=y*100l*100l+m*100l+d
                time=(h+mn/60.)*60.*60.
                k=where(data(*,3) eq date and $
                        data(*,2) gt time-60*60 and $
                        data(*,2) le time)
                if (k(0) ne -1) then begin
                    
                    for gi=0,n_elements(k)-1 do begin
                        h1=fix(data(k(gi),2)/(60.*60.))
                        mn1=fix((data(k(gi),2)-h1*60.*60.)/60.)
                        print, h1,mn1, data(k(gi),k_lon(0)),data(k(gi),k_pre(0))


                        if (data(k(gi),k_lat(0)) gt -180. and $
                            data(k(gi),k_lon(0)) gt -180. and $
                            data(k(gi),k_pre(0)) gt -180.) then begin
                            printf,1, counter,'146',d,m,y,h1,mn1,data(k(gi),k_lat(0)),$
                              data(k(gi),k_lon(0)), data(k(gi),k_pre(0)),$
                              format='(i5,x,a5,x,i2.2,"-",i2.2,"-",i4.4,x,i2.2,":",i2.2,x,f7.2,x,f7.2,x,f7.2 )'
                            counter=counter+1
                        endif
                    endfor
                endif
                for p=0,n_elements(sites)-1 do begin
                    
                    printf,1, counter,sites(p),d,m,y,h,mn,lats(p), lons(p), pres(p),$
                      format='(i5,x,a5,x,i2.2,"-",i2.2,"-",i4.4,x,i2.2,":",i2.2,x,f7.2,x,f7.2,x,f7.2 )'
                    counter=counter+1
                    
                endfor
            endfor
            
            
            printf,1,'99999   END  0- 0-   0  0: 0    0.00    0.00    0.00'
            
            close,1
        endfor
    endfor
endfor



end



