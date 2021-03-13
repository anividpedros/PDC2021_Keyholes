import rebound
import numpy as np
import matplotlib.pyplot as plt
import csv
import timeit
import re

sim = rebound.Simulation()

# Full Horizons data retrieval
eph = "JD2459184.5"
sim.add("Sun", date=eph)
sim.add("Mercury", date=eph)
sim.add("Venus", date=eph)
sim.add("Earth", date=eph) # Finds Earth-Moon Barycenter
sim.add("Mars", date=eph)
sim.add("Jupiter", date=eph)
sim.add("Saturn", date=eph)
sim.add("Uranus", date=eph)
sim.add("Neptune", date=eph)


# sim.add(id, m=0, date="JD2459184.5.0")
sim.add(m=0, a=1.,e=0.,inc=0.,Omega=0,omega=0,f=1.7501549, date=eph)
part = 9
earth = 3
id = "PDCasteroid"
#----------------------------------------------------
ps = sim.particles

Torb = 2.*np.pi # = 1 year
Noutputs = 1000
times = np.linspace(0, 1000.*Torb, Noutputs)

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
tplot = times/2/np.pi

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

