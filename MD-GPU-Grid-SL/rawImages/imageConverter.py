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

for j in range(200):
    i = j+1
    filename = "rawImage_s"+str(i)+".rawdata"
    
    data =  np.fromfile(filename,np.dtype('u1'),640*360*4)
    data = data.reshape(360,640,4);
    
    scipy.misc.toimage(data).save('conv'+str(i)+'.png')