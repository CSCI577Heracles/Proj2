import numpy as np
import math

DIST_CUTOFF = 2 ** (1 / 6.)


class Force(object):

    def __init__(self, c, GAMMA, GRAVITY):
        self.c = c
        self.GAMMA = GAMMA
        self.GRAVITY = GRAVITY

    def lj_force(self, mag, hat):
        eps = 1.0
        sig = 1.0
        a = (24 * eps) / mag * (2 * (sig / mag) ** 12) * hat
        a = np.nan_to_num(a)
        a = -a
        a[self.c.dr() > DIST_CUTOFF] = 0.
        return np.sum(a, axis=1)

    def damp_force(self, dvx, dvy, dx, dy, r_mag):
        gamma = self.GAMMA
        print "gamma: " + str(gamma)
        #return -gamma * (np.dot(v, r)) * (r / r ** 2)
        return gamma * (dvx * dx + dvy * dy) / r_mag

    def gravity(self):
        return np.ones(np.size(self.c.x)) * -self.GRAVITY * 0.01

    def a(self):
        dx = self.c.dx()
        dy = self.c.dy()
        #dz = self.c.dz()
        #dr = self.c.dr()

        dvx = self.c.dv_x()
        dvy = self.c.dv_y()

        #print "dvx"
        #print dvx[dvx > 0.]

        #print "dvy"
        #print dvy[dvy > 0.]

        r_mag = np.sqrt(dx ** 2 + dy ** 2)

        #print "r_mag"
        #print r_mag

        #r_mag = np.nan_to_num(r_mag)

        x_hat = dx / r_mag


        #x_hat = np.nan_to_num(x_hat)
        x_hat[-np.isfinite(x_hat)] = 0.

        #print "x_hat"
        #print x_hat

        y_hat = dy / r_mag
        #y_hat = np.nan_to_num(y_hat)
        y_hat[-np.isfinite(y_hat)] = 0.

        #print "y_hat"
        #print y_hat

        damp_force = self.damp_force(dvx, dvy, dx, dy, r_mag)
        #damp_force = np.nan_to_num(damp_force)
        damp_force[-np.isfinite(damp_force)] = 0.

        damp_force[self.c.dr() > DIST_CUTOFF] = 0.

        #print "damp force > 0"
        #print damp_force[damp_force > 0]
        #print "damp_force dimensions: " + str(damp_force.shape)
        ax = self.lj_force(r_mag, x_hat)
        ax += np.sum(np.nan_to_num(damp_force * x_hat), axis=1)


        # y accelerations:



        ay = self.lj_force(r_mag, y_hat)
        ay += np.sum(np.nan_to_num(damp_force * y_hat), axis=1)
        ay += self.gravity()

        ax[:self.c.NUM_SIDE * 2 + self.c.NUM_FLOOR] = 0.
        ay[:self.c.NUM_SIDE * 2 + self.c.NUM_FLOOR] = 0.

        #print "ax"
        #print ax

        #print "ay"
        #print ay

        return ax, ay


