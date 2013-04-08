import Container
import Force
import Integrator
import ContainerInitializer

import numpy as np
import matplotlib.pyplot as plt

NUM_TIMESTEPS = 4000
FRAME_RATE = 20
DELTA_T = 0.01
SQUEEZE = False
SQUEEZE_FACTOR = 0.997
NL = False                   # set True for Neighborlist, False for regular

if NL:
    print "Neighborlist ON"
else:
    print "Neighborlist OFF"

NL_UPDATE_RATE = 20
NL_DIST = 2. * 2. ** (1./6.)



state_list = []
pe_list = []
ke_list = []
pressure_list = []

def circle( xy, radius, color="lightsteelblue", facecolor="green", alpha=.6, ax=None ):

    e = plt.Circle( xy=xy, radius=radius )
    #if ax is None:
    ax = plt.gca()
    #ax.cla()
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_edgecolor( color )
    e.set_linewidth(3)
    e.set_facecolor( facecolor )
    e.set_alpha( alpha )

c = ContainerInitializer.ContainerInitializer("eight").getContainer()
c.init_nl()
c.update_nl(NL_DIST)

f = Force.Force(c, NL)
i = Integrator.Integrator(DELTA_T, f)

state_list = []
count = 0

plt.figure(1)
plt.clf()
plt.ion()
plt.xlim((0, 10))
plt.ylim((0, 10))
plt.grid()
ax = plt.gca()
plt.show()

while count < NUM_TIMESTEPS:
    #print "--------- BEGIN TIMESTEP " + str(count) + " --------------"

    if count > 3. and c.Lx > 8. * 2. ** (1/6.) and SQUEEZE:
        c.Lx *= SQUEEZE_FACTOR
        c.Ly *= SQUEEZE_FACTOR
        c.x *= SQUEEZE_FACTOR
        c.y *= SQUEEZE_FACTOR

    #if count % NL_UPDATE_RATE == 0:
    if True:
        c.update_nl(NL_DIST)

    i.integrate()
    print "--------------------------------------"
    print "Timestep: " + str(count*DELTA_T)
    print "aX:"
    print c.ax
    print "aY:"
    print c.ay
    pe_list.append(f.pe())
    ke_list.append(f.ke())

    #print "AX TIMESTEP " + str(count)
    #print c.ax

    #print "AY TIMESTEP " + str(count)
    #print c.ay

    #plt.plot(c.x, c.y)
    #plt.show()
    if count % FRAME_RATE == 0:
        #print "c.x"
        #print c.x
        #print "c.y"
        #print c.y
        plt.cla()
        for particle in range(c.x.size):
            part_x = 0.
            part_y = 0.
            circle((c.x[particle], c.y[particle]), radius = 0.5*2**(1/6.), ax=ax, facecolor='green')
        plt.draw()
        print "Timeunit: " + str(count * DELTA_T)
            #particles[particle] = {"x": c.x[particle], "y": c.y[particle]}
            #plt.xlim((0, 10))
            #plt.ylim((0, 10))
            #plt.plot(c.x[particle], c.y[particle], 'o')
        #plt.show()

    count += 1

time = np.linspace(0, NUM_TIMESTEPS*DELTA_T, NUM_TIMESTEPS)

plt.clf()
plt.plot(time, pe_list)
plt.ylabel('Potential Energy Per Particle')
plt.title('PE - Two Initial Condition')
plt.xlabel('Time Units')
#plt.savefig('eight_pe_nnl.png')
plt.show(block=True)

plt.clf()
plt.plot(time, pe_list)
plt.ylabel('Kinetic Energy')
plt.title('Kinetic Energy Per Particle For Shrinking Box')
plt.xlabel('Time Units')
#plt.savefig('prob2_pe.png')
plt.show(block=True)
