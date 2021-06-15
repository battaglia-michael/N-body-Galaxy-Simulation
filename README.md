# N-body-Galaxy-Simulation
Simulate an N-body galaxy using a Barnes-Hut recursive tree algorithm.
The speed for the simulation was then compared and contrasted against a particle mesh method and a brute force method O(N^2).
There is a full writeup in N-bodyReport.pdf, as well as a folder showing a few cool simulation demos. Below there are some explanations for each Python file extracted from the report. Was coded in Python 2.7 using numpy, scipy, matplotlib packages.


### BHSim.py
is a script that runs the Barnes-Hut simulation. It has various user controllable parameters such as number of particles, radius of the simulation domain, mass of the galaxy, scale length of the galaxy (density drops to 1/e at scale length), Barnes-Hut resolution (θBH), softening length, time step, and total simulation time. The script produces an animation showing time evolution of the system.

### Body.py 
defines the Body class, which is an object representing a physical particle in 2D space. It carries the particle’s mass (float, units of GMs), position (2D array, units of kpc), velocity (2D array, units of kpc/10Myr), force (2D array, units of GMs∗kpc/(10Myr)^2) box boundary (in kpc, for periodic boundary condition), and color (for plotting). Class carries functions for performing leapfrog steps (implementing equations in Equation 3.7), functions for computing distance between two body objects and checking position w.r.t. quadrants,
functions for computing kinetic energy and potential energy, mutators for force on body, and a function for plotting.

### Quad.py
defines the Quad class, which is an object representing a quadrant in 2D space. It carries the position of its lower left corner (2D array), and its side length (float). Class carries functions that return Quad objects representing each of its subquadrants, and a function for plotting.

### BHTree.py 
defines the BHTree class, which is an object representing a node of the Barnes-Hut tree. It carries a Quad object representing the quadrant that the node occupies. If it is not a stud node, it carries a Body object representing the aggregate particle in the node. If it has been given more than one particle, it carries 4 BHTree objects that are its child nodes. Class carries a function for inserting Bodies into the node. If it is a stub node, it accepts the Body and becomes a leaf node. If it is a leaf node, it changes its Body object into a aggregate particle, and gives its two Bodies to the appropriate children. If it is a branch node, it updates its aggregate Body and gives the new Body to its appropriate child. Class also carries a function for computing force of the node on a Body. If Body and the node’s aggregate Body satisfy Equation 3.1, then force is computed between the two. If not, then force is computed from each of the node’s children on the Body instead. If the node is a stub, then no force is added to the Body. Class also carries a function for recursively plotting the entire BHTree, with all non-stub children, grandchildren, etc

### MCgalaxy.py 
contains functions for generating various initial distributions
of particles using Monte Carlo methods. It can generate a galaxy using methods described in Section 2, and it can also generate a uniform distribution of particles with random velocities.

### BHtest.py 
is a script for testing speed and energy conservation of BarnesHut simulation.

### BruteForce.py
is a script for testing speed of O(N^2) brute force simulation
