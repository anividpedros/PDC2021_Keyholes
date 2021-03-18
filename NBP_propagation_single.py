import rebound
import numpy as np
import matplotlib.pyplot as plt
import csv
import timeit
import spiceypy as spice
import re


sim = rebound.Simulation()

# Full Horizons data retrieval
et = 719149020.237186 #epoch time

spice.furnsh('SPICEfiles/naif0012.tls.pc')
spice.furnsh('SPICEfiles/gm_de431.tpc')
spice.furnsh('SPICEfiles/pck00010.tpc')
spice.furnsh('SPICEfiles/2021_PDC-s11-merged-DE431.bsp')
# Change PATCH if other user
spice.furnsh('/Users/anivid/ExampleMICE/kernels/spk/de431_part-1.bsp')
spice.furnsh('/Users/anivid/ExampleMICE/kernels/spk/de431_part-2.bsp')

GM_sun = spice.bodvrd('SUN','GM',1)[1]
GM_mer = spice.bodvrd('199','GM',1)[1]
GM_ven = spice.bodvrd('299','GM',1)[1]
GM_ear = spice.bodvrd('399','GM',1)[1]
GM_mar = spice.bodvrd('4','GM',1)[1]
GM_jup = spice.bodvrd('5','GM',1)[1]
GM_sat = spice.bodvrd('6','GM',1)[1]
GM_ura = spice.bodvrd('7','GM',1)[1]
GM_nep = spice.bodvrd('8','GM',1)[1]

#SUN
state1,ls = spice.spkezr('199', et, 'ECLIPJ2000', 'None', 'SUN')
state1 = np.array(state1)
state2,ls = spice.spkezr('299', et, 'ECLIPJ2000', 'None', 'SUN')
state2 = np.array(state2)
state3,ls = spice.spkezr('399', et, 'ECLIPJ2000', 'None', 'SUN')
state3 = np.array(state3)
state4,ls = spice.spkezr('4', et, 'ECLIPJ2000', 'None', 'SUN')
state4 = np.array(state4)
state5,ls = spice.spkezr('5', et, 'ECLIPJ2000', 'None', 'SUN')
state5 = np.array(state5)
state6,ls = spice.spkezr('6', et, 'ECLIPJ2000', 'None', 'SUN')
state6 = np.array(state6)
state7,ls = spice.spkezr('7', et, 'ECLIPJ2000', 'None', 'SUN')
state7 = np.array(state7)
state8,ls = spice.spkezr('8', et, 'ECLIPJ2000', 'None', 'SUN')
state8 = np.array(state8)

# kep_earth = cspice_oscelt(state3,et,GM_sun)
sim.add(m=GM_sun)
sim.add(m=GM_mer, x=state1[0],y=state1[1],z=state1[2],vx=state1[3],vy=state1[4],vz=state1[5]) #Mercury
sim.add(m=GM_ven, x=state2[0],y=state2[1],z=state2[2],vx=state2[3],vy=state2[4],vz=state2[5])  #Venus
sim.add(m=GM_ear, x=state3[0],y=state3[1],z=state3[2],vx=state3[3],vy=state3[4],vz=state3[5]) #Earth
sim.add(m=GM_mar, x=state4[0],y=state4[1],z=state4[2],vx=state4[3],vy=state4[4],vz=state4[5]) #Mars
sim.add(m=GM_jup, x=state5[0],y=state5[1],z=state5[2],vx=state5[3],vy=state5[4],vz=state5[5]) #Jupiter
sim.add(m=GM_sat, x=state6[0],y=state6[1],z=state6[2],vx=state6[3],vy=state6[4],vz=state6[5]) #Saturn
sim.add(m=GM_ura, x=state7[0],y=state7[1],z=state7[2],vx=state7[3],vy=state7[4],vz=state7[5]) #Uranus
sim.add(m=GM_nep, x=state8[0],y=state8[1],z=state8[2],vx=state8[3],vy=state8[4],vz=state8[5]) #Neptune
     
          
# Asteroid PDC 2021
sim.add(m=0, a=333531592.021882,e=0.57830743096097,inc=0.313259712927621,Omega=3.64891542737495,omega=2.56612524483412,M=1.06084661066728)
part = 9
earth = 3
id = "PDCasteroid"
#----------------------------------------------------
ps = sim.particles

Torb = 365.25*3600*24 #1 year
Noutputs = 20001
times = np.linspace(0, 20.*Torb, Noutputs)

sim.integrator = "ias15" # IAS15 is the default
sim.move_to_com()    


d = np.zeros(Noutputs)
av= np.zeros(Noutputs)
ev= np.zeros(Noutputs)
iv= np.zeros(Noutputs)
Ov= np.zeros(Noutputs)
wv= np.zeros(Noutputs)
Mv= np.zeros(Noutputs)

start = timeit.default_timer()
for i,time in enumerate(times):
    sim.integrate(time)
    xa = ps[part].x
    ya = ps[part].y
    za = ps[part].z
    xt = ps[earth].x
    yt = ps[earth].y
    zt = ps[earth].z

    # distance computation
    d [i] = np.sqrt((xa-xt)**2 + (ya-yt)**2 + (za-zt)**2)
    av[i] = ps[part].a
    ev[i] = ps[part].e
    iv[i] = ps[part].inc
    Ov[i] = ps[part].Omega
    wv[i] = ps[part].omega
    Mv[i] = ps[part].f
stop = timeit.default_timer()
print ("Runtime: ", stop - start, " seconds")


#-----------------------
tplot = times

with open("PropagationFiles/"+id,'w') as f:
    f = csv.writer(f, delimiter=',')
    f.writerow(tplot)
    f.writerow(d)
    f.writerow(av)
    f.writerow(ev)
    f.writerow(iv)
    f.writerow(Ov)
    f.writerow(wv)
    f.writerow(Mv)

