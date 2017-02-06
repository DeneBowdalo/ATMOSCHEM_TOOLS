 
 for j=2012L,2012L do begin
     for i=0L,0L do begin
     date=j*10000L+(i+101)
     bpch_sep, 'ctm.bpch', 'ctm.bpch.'+string(date,format='(i8)'), tau0=nymd2tau(date)
     endfor 
     endfor

end
