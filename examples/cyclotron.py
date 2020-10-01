import matplotlib.pyplot as plt # Plottting module
from matplotlib.animation import FuncAnimation # Class for animation of solution
import sys
sys.path.append('../src')
from classes import *

B0 = 1 # B-field magnitude  [T]
qe = 1.602e-19 # electron charge [C]
me = 9.11e-31 # electron mass [kg]

solver = ParticleSolver(q=-qe, m=me)

def E(x,y,z):
    Ex, Ey, Ez = 0., 0., 0.
    return Ex, Ey, Ez

def B(x,y,z):
    Bx, By, Bz = 0., 0., B0
    return Bx, By, Bz

vx, vy, vz = [0., 1e4, 0.]
vperp = np.sqrt(vx**2 + vy**2)

omega_c = abs(solver.q/solver.m * B0)

T_c = 2*np.pi/omega_c

r_L = vperp/omega_c

x,y,z = [r_L, 0., 0.]

X0 = np.array([[x,y,z,vx,vy,vz]])

Nperiods = 1
Nsteps_per_period = 50
time = np.linspace(0., T_c*Nperiods, Nperiods*Nsteps_per_period)

solution = solver.solve(time, X0, E, B)
x = solution[0,:,0]/r_L
y = solution[0,:,1]/r_L
z = solution[0,:,2] 

fig,ax = plt.subplots()
l, = ax.plot([],[], marker='o')
ax.set_xlim(1.25*x.min(), 1.25*x.max())
ax.set_ylim(1.25*y.min(), 1.25*y.max())
ax.set_box_aspect(1.)
ax.set_aspect(1., adjustable='datalim')

def update(i):
    l.set_data(x[i],y[i])

ani = FuncAnimation(fig, update, range(x.size))
plt.show()
