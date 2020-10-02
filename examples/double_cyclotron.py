import matplotlib.pyplot as plt # Plottting module
from matplotlib.animation import FuncAnimation # Class for animation of solution
import sys # module needed to load in ../src/classes.py
sys.path.append('../src')
from classes import * # get all classes from ../src/classes.py

B0 = 1 # B-field magnitude  [T]
qe = 1.602e-19 # electron charge [C]
me = 9.11e-31 # electron mass [kg]
mp = 1.67e-27 # proton mass [kg]

# Make instance of Particle solver class and tell it that you are looking at an electron
esolver = ParticleSolver(q=-qe, m=me)
psolver = ParticleSolver(q=qe, m=mp)

# Define electric field [V]
def E(x,y,z):
    Ex, Ey, Ez = 0., 0., 0.
    return Ex, Ey, Ez

# Define magnetic field [T]
def B(x,y,z):
    Bx, By, Bz = 0., 0., B0
    return Bx, By, Bz

# initial velocites [m/s]
vx, vy, vz = [0., 1e4, 0.]
vperp = np.sqrt(vx**2 + vy**2) # velocity perpendicular to B-field

# cyclotron frequency [rad/s]
omega_ce = abs(esolver.q/esolver.m * B0)
omega_cp = abs(psolver.q/psolver.m * B0)

# cyclotron period [s]
T_ce = 2*np.pi/omega_ce
T_cp = 2*np.pi/omega_cp

# Larmor radius [m]
r_Le = vperp/omega_ce
r_Lp = vperp/omega_cp

# initial positions [m]
x,y,z = [0., 0., 0.]

# initial condition array
X0 = np.array([[x,y,z,vx,vy,vz]])

# Times to solve for
Nperiods = 1
Nsteps_per_period = 100
etime = np.linspace(0., T_ce*Nperiods, Nperiods*Nsteps_per_period)
ptime = np.linspace(0., T_cp*Nperiods, Nperiods*Nsteps_per_period)

# Solve the problem using Boris-Bunemann algorithm
esolution = esolver.solve(etime, X0, E, B)
psolution = psolver.solve(ptime, X0, E, B)

# Pull out positions for convience and normalize by Larmor radius
xe = esolution[0,:,0]
ye = esolution[0,:,1]
ze = esolution[0,:,2]
vxe = esolution[0,:,3]
vye = esolution[0,:,4]

xp = psolution[0,:,0]
yp = psolution[0,:,1]
zp = psolution[0,:,2]
vxp = psolution[0,:,3]
vyp = psolution[0,:,4]

# Plot trajectories of both particles
plt.plot(xe, ye)
plt.plot(xp, yp)

# Add arrows showing direction at time step 25
plt.arrow(xe[25], ye[25], 0.3*r_Le*vxe[25]/vperp, 0.3*r_Le*vye[25]/vperp, width = 0.005*r_Le, color='k')
plt.arrow(xp[25], yp[25], 0.3*r_Lp*vxp[25]/vperp, 0.3*r_Lp*vyp[25]/vperp, width = 0.005*r_Lp, color='k')

# make aspect ration 1. and display
plt.axis('equal')
plt.show()
