# -*- coding: utf-8 -*-
"""
Trajetória do míssil balístico Nodong

@author: José Rodrigues
"""

import matplotlib.pyplot as plt
import numpy as np
from earth import *
from nodong import *

# fig = plt.figure()
# g1 = fig.add_subplot(321)
# g1.set_aspect('auto')
# g1.set_xlabel('x [m]')
# g1.set_ylabel('y [m]')
# g2 = fig.add_subplot(323)
# g2.set_aspect('auto')
# g2.set_xlabel('t [s]')
# g2.set_ylabel('v [km/s]')
# g3 = fig.add_subplot(324)
# g3.set_aspect('auto')
# g3.set_xlabel('t [s]')
# g3.set_ylabel('at/g0')
# g5 = fig.add_subplot(325)
# g5.set_aspect('auto')
# g5.set_xlabel('t [s]')
# g5.set_ylabel('alph [º]')
# g4 = fig.add_subplot(322)
# g4.set_aspect('auto')
# g4.set_xlabel('t [s]')
# g4.set_ylabel('h [km]')
# g6 = fig.add_subplot(326)
# g6.set_aspect('auto')
# g6.set_xlabel('t [s]')
# g6.set_ylabel('teta [º]')


'CONDIÇÕES DE LANÇAMENTO'
pstart = 10                #Tempo do lift-off [s]
pduration = 30             #Tempo do Pitchover [s]        
pitchrate = 1              #Variaçao de alph em [º/s]

assert pitchrate <= 2, 'O pitchrate máximo é 2 [º/s]'

'CONDIÇÕES INICIAIS'
x = 0
y = RT + 2 #Raio da Terra [m]
h = 0
hmin = 1
t = 0          
dt = 0.001 # [s]
vx = 0
vy = 0         
alph = 0
psi = 90 # [º]

def rho(h):
    rho,p,t = atmosphere(h)
    return rho

def Fd(rho,vx,vy): #falta a direção!!!
    Fd = -0.5*Ca*S*rho(h)*(vx**2 + vy**2)
    return Fd

def Fthx(psi): 
    if t<=tb:
        Fthx = Isp*R*g0*(np.cos(radians(psi))) # o argumento de cos é em rad
    else:
        Fthx = 0
    return Fthx
def Fthy(psi):
    if t<=tb:
        Fthy = Isp*R*g0*(np.sin(radians(psi))) # o argumento de sin é em rad
    else:
        Fthy = 0
    return Fthy

def Fg(x,y): #falta a direção!!!
   # x += 1
   y >= RT
   Fg = -G*MT/(x**2+y**2)**1.5   
   return Fg

def m(t):
    m = m0 #massa inicial [kg]
    if t <= tb:
        m += -R*t
    else:
        m = me + payload
    return m

"""
TRAJETÓRIA
"""
while t >= 0 and h >= 0:
    if t < pstart:                                                             #"""LIFT-OFF"""                               
        ay = (Fthy(psi) + Fg(x,y)*y + Fd(rho,vx,vy)*vy)/m(t)
        vy += ay*dt
        y += vy*dt
        h = y - RT    
        t += dt
        assert sqrt(0+ay**2) <= ACCELMAX,'Missil explodiu: ultrapassou a aceleração máxima suportada'

    elif t < pstart + 0.5*pduration:                                           #"""PITCHOVER1"""    
        ax = (Fthx(psi) + Fg(x,y)*x + Fd(rho,vx,vy)*vx)/m(t)
        vx += ax*dt
        x += vx*dt
        ay = (Fthy(psi) + Fg(x,y)*y + Fd(rho,vx,vy)*vy)/m(t)
        vy += ay*dt
        y += vy*dt
        h = y - RT
        t += dt
        alph = pitchrate*(t-pstart)
        teta = degrees(acos(vx/(sqrt(vx**2+vy**2))))
        psi = teta - alph
        
    elif t < pstart + pduration:                                               #"""PITCHOVER2"""
        ax = (Fthx(psi) + Fg(x,y)*x + Fd(rho,vx,vy)*vx)/m(t)
        vx += ax*dt
        x += vx*dt
        ay = (Fthy(psi) + Fg(x,y) + Fd(rho,vx,vy)*vy)/m(t)
        vy += ay*dt
        y += vy*dt
        h = y - RT
        t += dt
        alphmax = pitchrate*0.5*pduration
        assert alphmax <= ALPHAMAX, 'Missil explodiu: Excedido o máximo do ângulo de Ataque'
        alph = alphmax - pitchrate*(t-(pstart + 0.5*pduration))
        teta = degrees(acos(vx/(sqrt(vx**2+vy**2))))
        psi = teta - alph
        
        assert sqrt(vx**2+vy**2)>VPMIN*1000/3600, 'Missil explodiu: não atingiu a velocidade mínima necessária'
        assert sqrt(vx**2+vy**2)<VPMAX*1000/3600, 'Missil explodiu: ultrapassou a velocidade máxima suportada'
        assert sqrt(ax**2+ay**2) <= ACCELMAX,'Missil explodiu: ultrapassou a aceleração máxima suportada'
        
    elif t > pstart + pduration:  	                                           #"""GRAVITY TURN"""                               
        ax = (Fthx(psi) + Fg(x,y)*x + Fd(rho,vx,vy)*vx)/m(t)
        vx += ax*dt
        x += vx*dt
        ay = (Fthy(psi) + Fg(x,y)*y + Fd(rho,vx,vy)*vy)/m(t)
        vy += ay*dt
        y += vy*dt
        h = y - RT
        t += dt    
        teta = degrees(acos(vx/(sqrt(vx**2+vy**2))))
        psi = teta #alph = 0
        
        assert sqrt(ax**2+ay**2) <= ACCELMAX,'Missil explodiu: ultrapassou a aceleração máxima suportada'
        
    #g3.scatter(t,sqrt(ax**2+ay**2)/g0,marker='o',c='blue',s=2.0)  
    #g1.scatter(x,y,marker='o',c='blue',s=2.0)         
    #g4.scatter(t,h,marker='o',c='red',s=2.0)
    print(Fthx(psi))
    assert h <= HMAX,'Missil explodiu: ultrapassou a altura máxima suportada' 
    
fi = degrees(acos(x/sqrt(x**2+y**2)))
dfi = 90 - fi
    
    #g3.scatter(t,sqrt(ax**2+ay**2)/g0,marker='o',c='blue',s=2.0)  
    #g1.scatter(x,y,marker='o',c='blue',s=2.0)         
    #g4.scatter(t,h,marker='o',c='red',s=2.0)

print('MISSION ACCOMPLISHED')
print('Valores de voo do missil Nodong:','\n','Alcance: ',RT*dfi/1000,'km')
print('Tempo de voo: ',t,'s') #'\n','Altura máxima: ',hmax)
print('Ângulo de ataque máximo: ',alphmax,'º')
# print('Velocidade máxima: ',vmax,'\n','Aceleração máxima: ',amax)
