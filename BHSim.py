#################################################################
# Name:     BHsim.py                                            #
# Authors:  Michael Battaglia                                   #
# Course:   Phy407                                              #
# Inst.:    Paul Kushner                                        #
# Date:     December 17, 2016                                   #
# Function: Program implements quadtree algorithm for computing #
#           forces on N bodies interacting under gravity.       #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#essential imports
from Body import Body
from Quad import Quad
from BHTree import BHTree
from MCgalaxy import generateGalaxy

#function: main
if __name__ == '__main__':
    #Milky Way parameters (default)
    r0 = 3 #kpc, scale length of galaxy
    m0 = 50.0 #10^9 solar mass, mass of galaxy
    #simulation space
    N = 1000 #number of particles
    L = 15.0 #half length of box, kpc
    #Barnes-Hut simulation resolution
    theta = 1.0
    epsilon = theta*L/np.sqrt(N)
    #time evolution parameters
    dt = 0.1 #10Myr
    T = 10.0 #10Myr
    steps = int(T/dt)
    #generate 1000 masses in 15kpc box
    bodies = generateGalaxy(r0, m0, N, L)
    #generate Barnes-Hut tree on original grid
    tree = BHTree(Quad(-L,-L,2*L))
    #populate tree with bodies from list
    for body in bodies:
        tree.insertBody(body)
    #calculate force on every body from tree and evolve leapfrog step
    for body in bodies:
        body.resetForce()
        tree.applyForce(body, theta, epsilon)
        #take a half time step
        body.leapFrog(dt)

    #make list of objects for plotting
    images = []
    #plotting setup
    fig = plt.figure()
    ax = plt.axes(xlim=(-L, L), ylim=(-L, L))
    #evolve N-body in time
    for i in range(steps):
        #computation counter
        print "Computing time step "+str(i+1)+"/"+str(steps)
        #generate Barnes-Hut tree on original grid
        tree = BHTree(Quad(-L,-L,2*L))
        #populate tree with bodies from list
        for body in bodies:
            tree.insertBody(body)
        #calculate force on every body from tree and evolve
        for body in bodies:
            body.resetForce()
            tree.applyForce(body, theta, epsilon)
            #take a time step
            body.update(dt)
        #append to list of objects for plotting
        position = np.array([body.r for body in bodies]).T
        scatter, = ax.plot(position[0], position[1], 'k.')
        images.append((scatter,))
                
    anim = animation.ArtistAnimation(fig, images, interval=100, blit=True)
    anim.save('BH-Nbody'+str(N)+'.mp4')
    plt.show()
