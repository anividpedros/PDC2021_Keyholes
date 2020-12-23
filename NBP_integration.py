import rebound
import numpy as np
import matplotlib.pyplot as plt
import csv
import timeit
import re
from multiprocessing import Pool


bodylist = []

with open('PHAsList.csv',newline='') as csvfile:
    bodies = list(csv.reader(csvfile,delimiter ='\n'))
    
def PropAsteroid(ab):
    a = ab[0]
    b = ab[-1]+1

    sim = rebound.Simulation()

    # Full Horizons data retrieval
    # "JD2448257.5"
    sim.add("Sun", date="JD2459184.5")
    sim.add("Mercury", date="JD2459184.5")
    sim.add("Venus", date="JD2459184.5")
    sim.add("Earth", date="JD2459184.5") # Finds Earth-Moon Barycenter
    sim.add("Mars", date="JD2459184.5")
    sim.add("Jupiter", date="JD2459184.5")
    sim.add("Saturn", date="JD2459184.5")
    sim.add("Uranus", date="JD2459184.5")
    sim.add("Neptune", date="JD2459184.5")
    
    for body in bodies[:][a:b]:

        id = body[0].split()[0]

        if re.search('[A-Z]', id):
            id = id[0:4] + ' ' + id[4:]

        # sim.add(id, m=0, date="JD2459184.5.0")
        sim.add(id, m=0, date="JD2459184.5")
        part = 9

        #----------------------------------------------------
        ps = sim.particles

        Torb = 2.*np.pi # = 1 year
        Noutputs = 10000
        times = np.linspace(0, 1000.*Torb, Noutputs)

        sim.integrator = "ias15" # IAS15 is the default
        sim.move_to_com()    

        x = np.zeros(Noutputs)
        y = np.zeros(Noutputs)
        av= np.zeros(Noutputs)
        ev= np.zeros(Noutputs)
        iv= np.zeros(Noutputs)
        Ov= np.zeros(Noutputs)
        wv= np.zeros(Noutputs)
        Mv= np.zeros(Noutputs)

        start = timeit.default_timer()
        for i,time in enumerate(times):
            sim.integrate(time)
            x[i] = ps[part].x
            y[i] = ps[part].y
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

        with open("Files/"+id,'w') as f:
            f = csv.writer(f, delimiter=',')
            f.writerow(tplot)
            f.writerow(av)
            f.writerow(ev)
            f.writerow(iv)
            f.writerow(Ov)
            f.writerow(wv)
            f.writerow(Mv)

        sim.remove(9)

cpus = 8
N = 2146
totalrange = range(0, N) # N # of simulations - ideally N multiple of cpus
ranges = np.array_split(totalrange, cpus)
pool = Pool(cpus)
pool.map(PropAsteroid,ranges)
