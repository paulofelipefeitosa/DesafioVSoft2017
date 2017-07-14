# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""

import numpy as np
import matplotlib.pyplot as plt
import math

i = 0
plt.figure(figsize=(12,12))
plt.ylim([0,500])
plt.xlim([0,500])
f = open('/home/paulofelipe/Códigos/Vsoft/VsoftSamplesDatabase/101_3.txt')
lin = int(f.readline())
for i in range(lin):
    x,y,orientation,type_pt = map(float,f.readline().split())
    plt.plot(x,y,'o')
    txt = 'P'+str(i)
    plt.text(x,y,txt)
    plt.arrow(x,y,math.cos(orientation)*15,math.sin(orientation)*15)
    i+=1
plt.savefig('three.png')

 
