# -*- coding: utf8 -*-
"""
Created on Tue Mar 21 03:30:36 2017

@author: Katsuhiro
"""

## rawImage To Image

import numpy as np
from matplotlib import pyplot as plt
import scipy.misc
#from Image 

width = 640*2
height = 360*2
wid_hei = width*height

for j in range(300):
    i = j+1
    filename = "rawImage_s"+str(i)+".rawdata"
    
    data =  np.fromfile(filename,np.dtype('u1'),wid_hei*4)
    data = data.reshape(height,width,4);
    
    scipy.misc.toimage(data).save('conv'+str(i)+'.png')