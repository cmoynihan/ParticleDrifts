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
    Ex, Ey, Ez = 0., 0., -100
    return Ex, Ey, Ez

# Define magnetic field [T]
def B(x,y,z):
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan(y/x)

    Bmag = B0/z
    # z = -Bmag(r-R0)^2
    Br =
    Bx, By, Bz = 0., 0., B0
    return Bx, By, Bz
