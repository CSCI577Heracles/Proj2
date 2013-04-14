import numpy as np
import networkx as nx 	# graph library - not needed but I think there is a fancy way to compute neighbor edges here

class Container(object):

    def __init__(self, DIST_CUTOFF):
        self._x = np.array([])
        self._y = np.array([])
        self._z = np.array([])

        self.vx = np.array([])
        self.vy = np.array([])
        self.vz = np.array([])

        self.ax = np.array([])
        self.ay = np.array([])
        self.az = np.array([])

        self.Lx = 0.
        self.Ly = 0.
        self.Lz = 0.

        self.DIST_CUTOFF = DIST_CUTOFF

    def init_nl(self):
        self.G = nx.DiGraph()
        for element in list(enumerate(self.x)):
            self.G.add_node(element[0])

    def update_nl(self, dist):
        self.G.remove_edges_from(self.G.edges())
        #print "dr:"
        #print self.dr()
        for row in list(enumerate(self.dr())):
            for col in list(enumerate(row[1])):
                if col[1] > dist:
                    self.G.add_edge(col[0], row[0], x=self.x[col[0]]-self.x[row[0]], y=self.y[col[0]]-self.y[row[0]], z=0., r=col[1])
        #print nx.adjacency_matrix(self.G)

    def neighbor(self, p):
        return self.G.neighbors(p)

    #def print_dy(self):
    #    for

    def x_dist(self, i, j):
        return self.G.edge[i][j]['x']

    def y_dist(self, i, j):
        return self.G.edge[i][j]['y']

    def r_dist(self, i, j):
        return self.G.edge[i][j]['r']

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x_prime):
        self._x = x_prime % self.Lx

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y_prime):
        self._y = y_prime % self.Ly

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, z_prime):
        #self._z = z_prime % self.Lz
        self._z = z_prime

    def add_particle(self, x, y, z, vx, vy, vz):
        self.x = np.hstack((self.x, x))
        #print "x: "
        #print self.x
        self.y = np.hstack((self.y, y))
        self.z = np.hstack((self.z, z))

        self.vx = np.hstack((self.vx, vx))
        self.vy = np.hstack((self.vy, vy))
        self.vz = np.hstack((self.vz, vz))
        #self.vx = vx
        #self.vy = vy
        #self.vz = vz

    def dx(self):
        #print "x: "
        #print self.x
        xtemp = np.tile(self.x, (self.x.size,1))
        #print "xtemp"
        #print xtemp
        dx = xtemp - xtemp.T
        dx[dx > self.Lx / 2.] -= self.Lx
        dx[dx < -self.Lx / 2.] += self.Lx

        #print "dx"
        #print dx
        return dx

    def dv_x(self):
        xtemp = np.tile(self.vx, (self.vx.size, 1))
        dvx = xtemp - xtemp.T

        return dvx

    def dv_y(self):
        ytemp = np.tile(self.vy, (self.vy.size, 1))
        dvy = ytemp - ytemp.T
        return dvy

    def dy(self):
        ytemp = np.tile(self.y, (self.y.size, 1))
        dy = ytemp - ytemp.T
        dy[dy > self.Ly / 2.] -= self.Ly
        dy[dy < -self.Ly / 2.] += self.Ly
        dy[dy > self.DIST_CUTOFF] = 0.
        return dy

    def dz(self):
        ztemp = np.tile(self.z, (self.z.size, 1))
        dz = ztemp - ztemp.T
        dz[dz > self.Lz / 2.] -= self.Lz
        dz[dz < -self.Lz / 2.] += self.Lz
        return dz

    def dr(self):
        # TODO: fix this, no negative values in dr matrix
        #print "dr:"
        #print np.sqrt(self.dx ** 2 + self.dy ** 2 + self.dz **2)
        return np.sqrt(self.dx() ** 2 + self.dy() ** 2)



    def dr2(self):
        r_mag = self.dx() ** 2 + self.dy() ** 2 + self.dz() ** 2
        return np.nan_to_num(r_mag)
