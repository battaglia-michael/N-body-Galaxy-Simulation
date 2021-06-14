#################################################################
# Name:     MCgalaxy.py                                         #
# Authors:  Michael Battaglia                                   #
# Course:   Phy407                                              #
# Inst.:    Paul Kushner                                        #
# Date:     December 17, 2016                                   #
# Function: Program generates N randomly distributed bodies     #
#           according to galactic density distribution and      #
#           naive rotation curve (no dark matter). To be used   #
#           as initializer for N-body simulation.               #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt

#essential imports
from Body import Body

#function: initialize array of bodies in a galaxy
def generateGalaxy(r0, m0, N, L):
    #divide mass of galaxy among N masses
    m = m0/N
    #generate N bodies
    bodies = []
    for i in range(N):
        #position from normalized distribution p(R) = (1/r0)exp(-r/r0)
        r = -r0*np.log(1.0-np.random.rand())
        if r < L:
            theta = 2.0*np.pi*np.random.rand()
            rx = r*np.cos(theta)
            ry = r*np.sin(theta)
            #velocity from naive estimate v ~ sqrt(GMgalaxy/r)
            v = 4.738*np.exp(-r0/r)/np.sqrt(r)
            vx = -v*np.sin(theta)
            vy = v*np.cos(theta)
            #generate body
            bodies.append(Body(m, rx, ry, vx=vx, vy=vy, L=L))
    return bodies

#function: initialize array of bodies in uniform box
def generateUniform(rho, v0, N, L):
    #divide total among N masses
    m = rho*L**2/N
    #generate N bodies randomly distributed in box
    bodies = []
    for i in range(N):
        vx, vy = np.random.rand()*2*v0-v0, np.random.rand()*2*v0-v0
        rx, ry = np.random.rand()*2*L-L, np.random.rand()*2*L-L
        #generate body
        bodies.append(Body(m, rx, ry, vx=vx, vy=vy, L=L))
    return bodies
    
