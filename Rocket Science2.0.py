# -*- coding: utf-8 -*-
"""
Trajetória do míssil balístico Nodong

Na linha de comandos colocar + 0 para ver apenas prints principais
Colocar + + para ver prints principais e secundarios
A analise dos melhores paramentros e feita no final (estao comentados)

@author: José Rodrigues
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from earth import *
from nodong import *
from sys import argv
# argv = ['','0','0']
print(argv)


#************* DEFINIÇÃO DE FUNCOES ******************
def flush_print_main(*args):
        if argv[1] == '+':
        # import datetime
        # now = datetime.datetime.now()
        # timestamp = now.strftime('%H:%M:%S.%f')
        # print (timestamp,  *args, flush=True)
            print (*args, flush=True)
            return
        
def flush_print_sec (*args):
        if argv[2] == '+':
        # import datetime
        # now = datetime.datetime.now()
        # timestamp = now.strftime('%H:%M:%S.%f')
        # print (timestamp,  *args, flush=True)
            print (*args, flush=True)
            return

def ac_g(pos):
    """
    Aceleracao gravitica

    pos: vetor posição
    NB.: np.linalg calcula o modulo do vetor
    """
    return -G*MT/np.linalg.norm(pos)**3 * pos

def Fd(h,vel): #Drag force
    """
    Forca de resistencia do ar

    h: altitude
    vel: vetor velocidade 
    """
    rho = atmosphere(h)[0]
    return -0.5*Ca*S*rho*np.linalg.norm(vel) * vel

def Fth(psi,t):
    """
    Força de impulsão

    psi: teta - alpha
    t: instante de tempo
    tburn: tempo de queima
    """
    if t <= tburn: 
        return Isp*ejectRate*g0*np.array([np.cos(psi),np.sin(psi)])
    else: 
        return np.array([.0,.0])


def Rocket(pstart = 10., pduration = 30., pitchrate = 1., plot = 0):
    """
    Realiza a simulacao do trajeto do foguetao

    pstart - duracao do liftoff
    pduration - duracao do pitchover
    pitchrate - taxa de varicao do alpha no pitchrate
    plot - opcao para incluir os graficos (0 ou 'p')
    """
    # print(pstart,pduration,pitchrate)
    # # ============ CONDICOES DE LANÇAMENTO =========================
    # pstart = 10                #Tempo do lift-off [s]
    # pduration = 30             #Tempo do Pitchover [s]        
    # pitchrate = 1.             #Variaçao de alph em [º/s]
    # #===============================================================
    pitchrate = np.radians(pitchrate)

    assert pitchrate <= np.radians(MAXPITCHRATE), \
        f'\nO pitchrate máximo foi ultrapassado:\
        \nPITCHRATEMAX -> {MAXPITCHRATE} \npitchrate -> {np.degrees(pitchrate)}'

    # ----- CONDICOES INICIAIS --------------------------------------------
    hInit = 2.                          # ~altura do CM
    posInit = np.array([.0,hInit + RT]) # vetor de posição inicial do CM (y=R_Terra + ~altura do CM)
    velInit = np.array([.0,.0])         # vetor de velocidade inicial do CM
    alphInit = .0                       # angulo de ataque inicial
    tetaInit = np.pi/2                  # teta inicial [rad]
    psiInit = tetaInit - alphInit       # psi inicial[rad]   
    dt = 0.001                          # passo de integracao [s]
    mfuel = mfuelInit
    # ---------------------------------------------------------------------

    t_list = [.0]
    height_list = [hInit]
    x_list,y_list = [posInit[0]],[posInit[1]]
    v_list = [.0]
    at_list = [.0]
    teta_list = [tetaInit]
    alph_list = [alphInit]
    psi_list = [psiInit]


    #ooooooooooooo0000000000000ºººººººººººº SIMULACAO ºººººººººººº0000000000000ooooooooooooo
    height = hInit
    pos = posInit
    vel = velInit
    t = 0

    while height >= 0:

        # MASSA DO ROCKET AO LONGO DA TRAJETÓRIA
        massRocket = me + payload + mfuel
        if mfuel > 0:
            mfuel -= ejectRate*dt
        else:
            mfuel = 0

        #"####### LIFTOF #######"
        if 0 <= t < pstart:

            alph = alphInit
            teta = tetaInit
            teta2 = tetaInit
            if t>pstart-dt: flush_print_sec(f'\nt = {t:.3f}s:',"LIFTOF COMPLETED")
        
        #"######## PITCHOVER #########"
        elif t <= pstart+pduration/2:

            if VPMIN > np.linalg.norm(vel) > VPMAX:
                flush_print_main(f'\nA velocidade do foguete tem que estar entre [{VPMIN:.3f},{VPMAX:.3f}].\
                        \nvel = {np.linalg.norm(vel):.3f}')
                break
            if alph > np.radians(ALPHAMAX):
                flush_print_main(f'\n O foguete ultrapassou o alpha máximo:\
                        \nALPHMAX -> {ALPHAMAX} \nalpha -> {np.degrees(alph)}')
                break
            
            alph += pitchrate*dt
            teta = np.arcsin(vel[1]/np.linalg.norm(vel))
            # teta2 = np.arccos(vel[0]/np.linalg.norm(vel))

        elif t <= pstart + pduration:

            alph -= pitchrate*dt
            teta = np.arcsin(vel[1]/np.linalg.norm(vel))
            # teta2 = np.arccos(vel[0]/np.linalg.norm(vel))
            if t>pstart+pduration-dt: flush_print_sec(f't = {t:.3f}s:',"PICHOVER COMPLETED")
        
        #"######## GRAVITY TURN ##########"
        else:
            alph = alphInit
            teta = np.arcsin(vel[1]/np.linalg.norm(vel))
            # teta2 = np.arccos(vel[0]/np.linalg.norm(vel))
            if tburn < t < tburn+dt: 
                flush_print_sec(f't = {t/60:.2f} min: FUEL FINISHED')
            if -0.01 < vel[1] < 0: flush_print_sec(f't = {t/60:.2f} min: BEGINNING DESCENDING MANOUVER')
        
        psi = teta - alph
        accel = ac_g(pos)+(Fd(height,vel)+Fth(psi,t))/massRocket
        vel += accel*dt
        pos += vel*dt
        height = np.linalg.norm(pos) - RT
        t += dt

        if height > HMAX:
            flush_print_main(f'\nO Foguete ultrapassou a altura máxima:\
                    \nHMAX -> {HMAX} \naltura -> {height:.3f}')
            break
        if np.linalg.norm(accel) > ACCELMAX:
            flush_print_main(f'\nO Foguete ultrapassou a aceleração máxima:\
                    \nACCELMAX -> {ACCELMAX} \naccel -> {np.linalg.norm(accel):.3f}')
            break
        
        t_list.append(t)
        height_list.append(height)
        x_list.append(pos[0])
        y_list.append(pos[1])
        v_list.append(np.linalg.norm(vel))
        teta_list.append(teta)
        alph_list.append(alph)
        psi_list.append(psi)
    #ooooooooooooo0000000000000ºººººººººººººººººººººººº0000000000000ooooooooooooo

    if height < 0:
        at_list = np.gradient(v_list,dt) #cálculo da aceleracao tangencial

        #*********** PRINTS FINAIS *****************
        flush_print_main(f't = {t/60:.2f} min: NODONG HAS LANDED SUCCESSFULLY WITH A PLAYLOAD OF {payload} kg!')
        flush_print_main(f'\n-> Range: {x_list[-1]/1000:.2f} km \n-> Time of Flight: {t_list[-1]/60:.2f} min \n-> Apogee: {max(height_list):.2f} km')
        flush_print_main(f'-> Maximum speed: {max(v_list)/1000:.2f} km/s (mach {int(max(v_list)/343)})')
        flush_print_main(f'-> Maximum acceleration: {max(at_list)/g0:.2f} g\n')

        if plot == 'p':
            x_earth_list = np.linspace(-1000,1500,len(x_list))
            y_earth_list = np.sqrt((RT/1000)**2-np.array(x_earth_list)**2)

            #*********** GRAFICOS *****************
            altura = 718 - 66
            largura = 1420 - 70
            px = 1/plt.rcParams['figure.dpi']  # pixel in inches
            figsize_init = (largura*px,altura*px)

            fig = plt.figure(figsize=figsize_init)
            fig.suptitle(
                f'Nodong-1 Simulation: pstart = {pstart}s, pduration = {pduration}s, pitchrate = {np.degrees(pitchrate):.1f}'+ r'$\circ/s$ '
                , fontweight="bold")
            
            plt.rcParams['axes.grid'] = True
          
            gs = GridSpec(4, 2, figure=fig)

            ax = fig.add_subplot(gs[0, :])
            ax.plot(np.array(x_list)/1000,np.array(y_list)/1000,c='C1',label='Rocket trajectory')
            ax.plot(x_earth_list,y_earth_list,c='b',label='Earth',linestyle=':')
            ax.set_ylabel('y [km]')
            ax.set_xlabel('x [km]')
            ax.set_xlim(-1000,1500)
            ax.legend(loc='upper right')

            ax = fig.add_subplot(gs[1, :])
            ax.plot(np.array(t_list)/60,(np.array(height_list))/1000,c='C1')
            ax.set_ylabel('height [km]')
            ax.set_xlabel('time [min]')

            ax1 = fig.add_subplot(gs[2, 0])
            ax1.plot(np.array(t_list)/60,np.array(v_list)/1000,c='C0')
            plt.setp(ax1.get_xticklabels(), visible=False)
            ax1.set_ylabel('v [km/s]')

            ax2 = fig.add_subplot(gs[2, 1],sharex=ax1)
            ax2.plot(np.array(t_list)/60,np.array(at_list)/g0,c='C0')
            plt.setp(ax2.get_xticklabels(), visible=False)
            ax2.set_ylabel('$a_t/g_0$')

            ax3 = fig.add_subplot(gs[3, 0])
            ax3.plot(np.array(t_list)/60,np.degrees(alph_list),c='C2')
            ax3.set_ylabel(r'$\alpha ~ [\circ]$')
            ax3.set_xlabel('t [min]')

            ax4 = fig.add_subplot(gs[3, 1],sharex=ax3)
            ax4.plot(np.array(t_list)/60,np.degrees(teta_list),c='C2')
            ax4.set_ylabel(r'$\theta ~ [\circ]$')
            ax4.set_xlabel('t [min]')

            gs.tight_layout(fig, rect=[0, 0, 1, 0.97])
            # gs.tight_layout(fig)
            plt.show()

        return max(x_list)
    
    else:
        return 0



############################## CÁLCULO DOS MELHORES RESULTADOS ####################################
# pstart_list = [i for i in range(3,18,3)]
# pduration_list = [i for i in range(5,45,5)]
# pitchrate_list = [i/10 for i in range(5,25,5)]

# results = {} 
# for pitc in pitchrate_list:
#     for pdu in pduration_list:
#         for ps in pstart_list:
#             # print(f'\n{ps} {pdu} {pitc}')
#             simulation = Rocket(ps,pdu,pitc)/1000
#             # print(f'xmax = {simulation}')
#             results[(ps,pdu,pitc)] = simulation

# results = dict(sorted(results.items(), key=lambda item: item[1], reverse=True)) #ordenar pelo maior alcance
# print(results)
# with open('RocketScience_Results.txt', 'w') as file:
#     # Write the dictionary to the file
#     file.write(str(results))

# bestComb = next(iter(results)) #melhor (primeiro) resultado de results
# print('\nMELHOR COMBINACAO:')
# print(f'pstart = {bestComb[0]}, pduration = {bestComb[1]}, pitchdrate = {bestComb[2]}')

# Rocket(bestComb[0], bestComb[1], bestComb[2], 'p')
########################################################################################

# BEST
Rocket(3,15,0.5,'p')

# FUNCIONA
# Rocket(10., 30., 1, 'p')