InFile = 'ctm.bpch_1'
 for j=2006L,2006L do begin
     date=j*10000L+101L
     date2=(j+1)*10000L+101L
     bpch_sep, InFile, 'ctm.bpch.'+string(date,format='(i8)'), tau0=nymd2tau(date), tau1=nymd2tau(date2)
 endfor

end
                                                                                                                
                                                                                                                
         
