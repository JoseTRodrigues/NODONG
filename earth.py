"""This module implements data and functions related to our
planet, Earth
"""

from math import exp

RT = 6378165. # earth radius"
MT = 5.974e24 # earth mass
G = 6.673e-11 # g gravitational constant
g0 = G*MT/RT**2 # g at earth surface

PS = 101325. # standard pressure at sea level
TS = 288.16 # standard temperature at sea level
LE = 0.0065 # temperature decrease per meter height in troposphere

RG = 287.26 # R divided by molar mass of air
AS = g0/(LE*RG)

HT = 11000.0 # height where troposphere ends 
TT = TS -HT*LE # temperature at end of troposphere

AT = g0/(TT*RG)
PT = PS*(TT/TS)**AS #standard pressure at end of troposphere

def atmosphere(h):
    """
    "Returns density of air (kg/m^3),pressure (Pa) and temperature (K) as 
    function of height(h) from earth surface (m)."
    """
    if h < HT:
        temp = TS-LE*h
        p = PS*(temp/TS)**AS
        rho = p/(RG*temp)
    else:
        p = PT*exp(AT*(HT-h))
        rho = p/(RG*TT)
        temp = TS-LE*HT
    return rho,p,temp

def atmosphere_profile(hmax=200000,hstep=1000,filename="atmosphere.dat"):
    """
    "Prints a table of h,rho(h),p(h)
    """
    h = 0.

    with open(filename,'w') as f:
        while h <= hmax:
            rho, p, temp = atmosphere(h)
            print("{:10.3f} {:10.3g} {:10.3g} {:10.3f} ".format(h,rho,p,temp),file=f)
            h += hstep



if __name__ == "__main__":
    atmosphere_profile()
    


