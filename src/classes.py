import numpy as np

class ParticleSolver:

    def __init__(self, q=-1.60217662e-19, m=9.10938356e-31, freq_correct=True):
        self.freq_correct = freq_correct
        self.q = q
        self.m = m

    def frequency_correction(self, x):
        if x == 0.:
            return 1.
        return np.tan(x)/x

    def velocity_push_back(self, X0, dt, E, B):
        qmdt2 = self.q/self.m*dt/2

        x  = X0[0]
        y  = X0[1]
        z  = X0[2]
        vx = X0[3]
        vy = X0[4]
        vz = X0[5]

        Ex, Ey, Ez = E(x,y,z)
        Bx, By, Bz = B(x,y,z)

        if self.freq_correct:
            alphax = self.frequency_correction(qmdt2*Bx)
            alphay = self.frequency_correction(qmdt2*By)
            alphaz = self.frequency_correction(qmdt2*Bz)
        else:
            alphax, alphay, alphaz = 1.,1.,1.

        # Half Acceleration (E-field)
        vx += qmdt2 * Ex * alphax
        vy += qmdt2 * Ey * alphay
        vz += qmdt2 * Ez * alphaz
        # Rotation (B-field)
        tx = qmdt2*Bx * alphax
        ty = qmdt2*By * alphay
        tz = qmdt2*Bz * alphaz
        tsq = tx*tx + ty*ty + tz*tz
        sx = 2*tx/(1+tsq)
        sy = 2*ty/(1+tsq)
        sz = 2*tz/(1+tsq)
        vpx = vx + vy*tz-vz*ty
        vpy = vy + vz*tx-vx*tz
        vpz = vz + vx*ty-vy*tx
        vx += vpy*sz-vpz*sy
        vy += vpz*sx-vpx*sz
        vz += vpx*sy-vpy*sx
        # Half Acceleration (E-field)
        vx += qmdt2 * Ex * alphax
        vy += qmdt2 * Ey * alphay
        vz += qmdt2 * Ez * alphaz
        # Push Position
        x += vx*dt
        y += vy*dt
        z += vz*dt
        return vx,vy,vz

    def boris_bunemann(self, time, X0, E, B):
        dt = time[1]-time[0]
        qmdt2 = self.q/self.m*dt/2

        N = np.size(time)
        M = np.size(X0)
        X = np.empty((N,M))
        X[:] = np.nan

        X[0,:] = X0

        x = X[0,0]
        y = X[0,1]
        z = X[0,2]
        vx = X[0,3]
        vy = X[0,4]
        vz = X[0,5]
        for n in range(0, N-1):
            Ex, Ey, Ez = E(x,y,z)
            Bx, By, Bz = B(x,y,z)

            if self.freq_correct:
                alphax = self.frequency_correction(qmdt2*Bx)
                alphay = self.frequency_correction(qmdt2*By)
                alphaz = self.frequency_correction(qmdt2*Bz)
            else:
                alphax, alphay, alphaz = 1.,1.,1.

            # Half Acceleration (E-field)
            vx += qmdt2 * Ex * alphax
            vy += qmdt2 * Ey * alphay
            vz += qmdt2 * Ez * alphaz
            # Rotation (B-field)
            tx = qmdt2*Bx * alphax
            ty = qmdt2*By * alphay
            tz = qmdt2*Bz * alphaz
            tsq = tx*tx + ty*ty + tz*tz
            sx = 2*tx/(1+tsq)
            sy = 2*ty/(1+tsq)
            sz = 2*tz/(1+tsq)
            vpx = vx + vy*tz-vz*ty
            vpy = vy + vz*tx-vx*tz
            vpz = vz + vx*ty-vy*tx
            vx += vpy*sz-vpz*sy
            vy += vpz*sx-vpx*sz
            vz += vpx*sy-vpy*sx
            # Half Acceleration (E-field)
            vx += qmdt2 * Ex * alphax
            vy += qmdt2 * Ey * alphay
            vz += qmdt2 * Ez * alphaz
            # Push Position
            x += vx*dt
            y += vy*dt
            z += vz*dt
            # Store the coordinates in X
            X[n+1,0] = x
            X[n+1,1] = y
            X[n+1,2] = z
            X[n+1,3] = vx
            X[n+1,4] = vy
            X[n+1,5] = vz
        return X


    def solve(self, time, X0, E, B):
        dt = time[1]-time[0]
        solution = np.zeros((len(X0), len(time), 6))
        for i in range(len(X0)):
            vx0,vy0,vz0 = self.velocity_push_back(X0[i,:], -dt/2., E, B)
            X0p = np.array( [X0[i,0], X0[i,1], X0[i,2], vx0, vy0, vz0] )
            solution[i,:] = self.boris_bunemann(time, X0p, E, B)
        return solution
