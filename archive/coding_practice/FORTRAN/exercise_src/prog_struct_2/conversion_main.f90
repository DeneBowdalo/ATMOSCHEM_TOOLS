PROGRAM conversion_main
Use conversion

Implicit None
Real(wp) :: value, degs, rads, con_hours, con_minutes, con_seconds
character(1) :: type_con 

print*,'Type d to convert from radians to degrees'
print*,'Type r to convert from degrees to radians'
print*, 'Type t to convert seconds into, hours, minutes and seconds'
read(*,*) type_con

print*,'Please enter a value to convert'
read(*,*) value

if(type_con == 'd') then
    call rad_to_deg(value, degs)
 
else if(type_con == 'r') then
    call deg_to_rad(value, rads)

else if(type_con == 't') then
    call convert_time(value, con_hours, con_minutes, con_seconds) 

endif
END PROGRAM conversion_main
