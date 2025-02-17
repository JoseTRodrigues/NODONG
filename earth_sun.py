#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 10:35:27 2019
Simulation of Earth movement around the Sun.
Integration using the Euler-Cromer algorithm.
Uses astronomic units.
@author: jap
"""

from math import *
import matplotlib.pyplot as plt

G = 4*pi**2 # gravitation constant in au

dt = 1./(365*24*60)

day = 24*60

tmax = 1.0
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 10:35:27 2019
Simulation of Earth movement around the Sun.
Integration using the Euler-Cromer algorithm.
Uses astronomic units.
@author: jap
"""

import numpy as np
import matplotlib.pyplot as plt

G = 4*np.pi**2 # gravitation constant in au

dt = 1./(365*24*60)

day = 24*60

tmax = 1.0

# initial earth velocity in a.u.

vx = 0.0
vy = 2*np.pi
v = np.array([vx,vy])

x = 1.0
y = 0.0
r = np.array([x,y])

t = 0.0

n = 0

plt.title('Earth orbit')
plt.axis('equal')
plt.xlabel('x/a.u.')
plt.ylabel('y/a.u.')

plt.scatter(0,0,s=100,color='orange') #draw Sun, circle of size 100 pixels

while t <= tmax:
    # ax = -G*x/(x**2+y**2)**1.5
    # ay = -G*y/(x**2+y**2)**1.5
    a = -G*r/np.linalg.norm(r)**3
    
    # vx +=  ax*dt  #same as vx = vx + ax*dt   
    # vy +=  ay*dt
    v += a*dt
    
    # x +=  vx*dt
    # y +=  vy*dt
    r += v*dt
    
    t += dt 
    n += 1
    
    if (n%day==0): 
        # plt.scatter(x,y,s=5, color='blue')
        plt.scatter(r[0],r[1],s=5, color='blue')
        # plt.pause(0.05) # omit this line to get a fast plot of the trajectory (no animation)

plt.show()
