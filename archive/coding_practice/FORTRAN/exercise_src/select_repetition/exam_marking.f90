PROGRAM marking
IMPLICIT None

integer ::ave, candidate, first_mark, second_mark, third_mark

do 
write(*,*)'Input a Candidate Number, enter a -ve '
read(*,*) candidate
if(candidate <0) EXIT

write(*,*)'Input First Mark'
read(*,*) first_mark
write(*,*)'Input Second Mark'
read(*,*) second_mark
write(*,*)'Input Third Mark'
read(*,*) third_mark

!Calculate Average
ave = nint((first_mark+second_mark+third_mark)/3)

if(ave >= 90) then
    write(*,*) 'Candiate Number', candidate, 'acheived an A'

else if(ave<90.and.ave>=85) then
    write(*,*) 'Candiate Number', candidate, 'acheived an AB'

else if(ave<85.and.ave>=80) then
    write(*,*) 'Candiate Number', candidate, 'acheived a B'

else if(ave<80.and.ave>=75) then
    write(*,*) 'Candiate Number', candidate, 'acheived a BC' 

else if(ave<75.and.ave>=70) then
    write(*,*) 'Candiate Number', candidate, 'acheived a C' 

else if(ave<70.and.ave>=65) then
    write(*,*) 'Candiate Number', candidate, 'acheived a CD' 

else if(ave<65.and.ave>=60) then 
    write(*,*) 'Candiate Number', candidate, 'acheived a D' 

else if(ave < 60) then
    write(*,*) 'Candiate Number', candidate, 'acheived an F'

end if
end do

END PROGRAM marking 
