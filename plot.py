# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 08:49:26 2017

@author: jorge
"""

plt.figure(1,figsize=(20, 10))
plt.clf()

plt.subplot(211)
plt.title("Cputime=%.2fs (%.2f/TimeStep)" % ((time.time()-time0),(time.time()-time0)/kiter))
plt.plot(tiempos,zmaxdet, 'bo-')
plt.ylabel('Zmin')
plt.grid()
plt.subplot(212)
plt.plot(tiempos,alphamin, 'r.-')
plt.plot(tiempos,alphamax, 'b.-')
plt.ylabel('Alpha Min/Max')
plt.xlabel('Tiempo')
plt.grid()
plt.pause(0.005)
#plt.show()
