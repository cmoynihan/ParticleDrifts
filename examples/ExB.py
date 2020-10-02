import matplotlib.pyplot as plt # Plottting module
from matplotlib.animation import FuncAnimation # Class for animation of solution
import sys # module needed to load in ../src/classes.py
sys.path.append('../src')
from classes import * # get all classes from ../src/classes.py

B0 = 1 # B-field magnitude  [T]
qe = 1.602e-19 # electron charge [C]
me = 9.11e-31 # electron mass [kg]

# Make instance of Particle solver class and tell it that you are looking at an electron
solver = ParticleSolver(q=-qe, m=me)

# Define electric field [V]
def E(x,y,z):
    Ex, Ey, Ez = 0., 300., 0.
    return Ex, Ey, Ez

# Define magnetic field [T]
def B(x,y,z):
    Bx, By, Bz = 0., 0., B0
    return Bx, By, Bz

# initial velocites [m/s]
vx, vy, vz = [0., 1e4, 0.]
vperp = np.sqrt(vx**2 + vy**2) # velocity perpendicular to B-field

# cyclotron frequency [rad/s]
omega_c = abs(solver.q/solver.m * B0)

# cyclotron period [s]
T_c = 2*np.pi/omega_c

# Larmor radius [m]
r_L = vperp/omega_c

# initial positions [m]
x,y,z = [r_L, 0., 0.]

# initial condition array
X0 = np.array([[x,y,z,vx,vy,vz]])

# Times to solve for
Nperiods = 5
Nsteps_per_period = 50
time = np.linspace(0., T_c*Nperiods, Nperiods*Nsteps_per_period)

# Solve the problem using Boris-Bunemann algorithm
solution = solver.solve(time, X0, E, B)

# Pull out positions for convience and normalize by Larmor radius
x = solution[0,:,0]/r_L
y = solution[0,:,1]/r_L
z = solution[0,:,2]/r_L

# Create figure for animation of motion
fig,ax = plt.subplots()
N = 10
lines = []
for i in range(N):
    lines.append(ax.plot([],[], marker='o')[0])

l, = ax.plot([],[], '-') # empty line that will be modified at each step
s, = ax.plot([],[], marker='o')
ax.set_ylim(1.25*y.min(), 1.25*y.max()) # set y limits
ax.set_xlim(1.25*x.min(), 1.25*x.max()) # set x limits
ax.axvline(-1) # vertical line at -r_L
ax.set_xlabel('$x/r_L$')
ax.set_ylabel('$y/r_L$')
ax.set_title(r'$E \times B$ Motion w\ $E = E_0\hat{y}$, $B = B_0 \hat{z}$')

# make aspect ratio 1. so circles appear correctly
ax.set_aspect(1., adjustable='datalim')

# function called to update plot at each time step
# this just sets the point to the location of the particle
def update(i):
    l.set_data(x[:i+1],y[:i+1])
    s.set_data(x[i],y[i])

# make the animation and display
ani = FuncAnimation(fig, update, range(x.size), interval = 50)
plt.show()
