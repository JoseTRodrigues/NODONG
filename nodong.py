#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 14:32:14 2020
Nodong-1 Missile data
@author: jap
"""
from math import log,radians,pi
from earth import g0

Isp = 226. #specific impulse (s)
mfuelInit = 12912. # fuel mass (kg)
tburn = 110. # burn out time (s), also known as t1
ejectRate = mfuelInit/tburn # rate of fuel burn (kg/s)
payload = 500 # payload (weapons) mass (kg)
me = 1780 #empty mass, no fuel nor payload
m0 = me + mfuelInit + payload # total starting mass
d = 1.32 # diameter of rocket (m)
S = pi*(d/2)**2 # cross-section of rocket (m^2)
Ca = 0.25 # aerodinamic coeficient for airdrag formula (dimensionless)
massr = (m0-mfuelInit)/m0 # mass ratio at take off
deltav = Isp*g0*log(1./massr) # expected delta-v (m/s, free space) from rocket equation

#Nodong physical limits
ALPHAMAX = 20 # maximum attack angle
VPMIN = 50.0  # minimum velocity for pitchover
VPMAX = 800.0 # maximum velocity for pitchover
MAXPITCHRATE = 2.0  #  maximum pitchover rate [º]
HMAX = 250.0E3  # maximum height for missile
ACCELMAX = 15*9.8 # maximum accelration for missile

