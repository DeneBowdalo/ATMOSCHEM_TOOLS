import numpy as np
import matplotlib.pyplot as plt

#acts as area
XA1 = -115.
XA2 = -80.
YA1 = -4.
YA2 = 3.

#acts as gridbox
XB1 = -110.
XB2 = -108.
YB1 = -6.
YB2 = 5.

Aheight = YA2-YA1
Bheight = YB2-YB1
Awidth = XA2-XA1
Bwidth = XB2-XB1

#how much of gridbox is covered by area
#B all in A
if (XB1 >= XA1) & (XB2 <= XA2):
    overlapfracx = 1. 
#PARTIAL B IN A
elif (XB1 >= XA1) & (XB1 < XA2) & (XB2 >= XA2):
    overlapx = XA2 - XB1    
    overlapfracx = overlapx/Bwidth
elif (XB1 < XA1) & (XB2 > XA1):
    overlapx = XB2 - XA1    
    overlapfracx = overlapx/Bwidth
else:
    overlapfracx = 0. 

#B all in A
if (YB1 >= YA1) & (YB2 <= YA2):
    overlapfracy = 1. 
#PARTIAL B IN A
elif (YB1 >= YA1) & (YB1 < YA2) & (YB2 >= YA2):
    overlapy = YA2 - YB1    
    overlapfracy = overlapy/Bheight

elif (YB1 < YA1) & (YB2 > YA1):
    overlapy = YB2 - YA1    
    overlapfracy = overlapy/Bheight
else:
    overlapfracy = 0.

overlap_ratio = overlapfracx*overlapfracy

plt.plot([XA1],[YA1],marker='x')
plt.plot([XA2],[YA2],marker='x')
plt.plot([XB1],[YB1],marker='o')
plt.plot([XB2],[YB2],marker='o')

plt.xlim(-5,21)
plt.ylim(-5,21)

plt.show()
