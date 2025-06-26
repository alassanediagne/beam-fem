import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.linalg as la


# Hermite basis polynomials on [-1,1]
def N1(x: float) -> float:
    return 0.5 - 0.75 * x + +0.25 * x**3


def N2(x: float) -> float:
    return 1 / 4 - 1 / 4 * x - 1 / 4 * x**2 + 1 / 4 * x**3


def N3(x: float) -> float:
    return 1 / 2 + 3 / 4 * x - 1 / 4 * x**3


def N4(x: float) -> float:
    return -1 / 4 - 1 / 4 * x + 1 / 4 * x**2 + 1 / 4 * x**3


# Hermite basis polynomials after change of variables onto [a,b]
def Nk1(x, a, b):
    xi = 2 / (b - a) * (x - a) - 1
    return N1(xi)


def Nk2(x, a, b):
    xi = 2 / (b - a) * (x - a) - 1
    return (b - a) / 2 * N2(xi)


def Nk3(x, a, b):
    xi = 2 / (b - a) * (x - a) - 1
    return N3(xi)


def Nk4(x, a, b):
    xi = 2 / (b - a) * (x - a) - 1
    return (b - a) / 2 * N4(xi)


def Kloc(L: float, E: float, I: float) -> np.ndarray:
    """
    returns stiffness matrix on one element
    """
    K = np.array(
        [
            [12, 6 * L, -12, 6 * L],
            [6 * L, 4 * L**2, -6 * L, 2 * L**2],
            [-12, -6 * L, 12, -6 * L],
            [6 * L, 2 * L**2, -6 * L, 4 * L**2],
        ]
    )
    K = (E * I / L**3) * K
    return K


def zero_bc(K, F, js):
    """
    applies zero boundary conditions to the system
    K: stiffness matrix
    F: left hand side vector
    j: indeces of the solution array that are set to zero
    returns reduced stiffness matrix and second member
    """
    K_red = K
    F_red = F
    K_red = np.delete(K_red, js, axis=0)
    K_red = np.delete(K_red, js, axis=1)
    F_red = np.delete(F_red, js)
    return K_red, F_red


def lower_banded(K):
    """
    transforms K to a lower banded matrix
    """
    m, _ = K.shape
    K_banded = np.zeros((4, m))
    K_banded[0, :] = np.diag(K)
    K_banded[1, :-1] = np.diag(K, 1)
    K_banded[2, :-2] = np.diag(K, 2)
    K_banded[3, :-3] = np.diag(K, 3)
    return K_banded


class DiscretizedBeam:
    def __init__(self, nodes, E, I):
        """
        Initialize beam with list of nodes, Young's modulus and moment of intertia
        """
        self.nodes: np.ndarray = np.array(nodes)
        self.nnodes: int = len(nodes)
        self.E: float = E
        self.I: float = I
        self.loads: np.ndarray = np.zeros(
            (2, len(nodes))
        )  # array with the nodal vertical point loads in the first row and nodal moment loads in the second row
        self.free: np.ndarray = np.ones(
            (2, len(nodes)), dtype=bool
        )  # array indicating whether the element is free or fixed

    def add_support(self, coordinate, type):
        """
        add support at node given by x-coordinate

        :param coordinate: x-coordinate of support
        :param type: is either 'fixed', 'Roller', 'hinge' or 'simple'
        """
        node = np.searchsorted(self.nodes, coordinate)
        if type == "fixed":
            self.free[0, node] = False
            self.free[1, node] = False
        elif type == "Roller":
            self.free[0, node] = False
        elif type == "Hinge":
            self.free[0, node] = False
        elif type == "simple":
            self.free[0, node] = False
        else:
            raise TypeError("Please enter an adequate support")

    def find_element(self, x):
        diff = np.subtract(x, self.nodes)
        k = np.argmin(np.where(diff >= 0, diff, np.inf))
        if x == self.nodes[-1]:
            k = self.nnodes - 2
        return k

    def add_pointload(self, coordinate, load):
        """
        :param (float) coordinate: x-coordinate where load is applied. Does not have to be a node
        :param (float) load: magnitude of load

        adds a concentrated point load.
        """
        if coordinate in self.nodes:
            self.loads[0,]
        k = self.find_element(coordinate)
        L_elt = self.nodes[k + 1] - self.nodes[k]
        a = coordinate - self.nodes[k]
        b = L_elt - a
        self.loads[0, k] -= load * b**2 * (L_elt + 2 * a) / L_elt**3
        self.loads[1, k] -= load * a * b**2 / L_elt**2
        self.loads[0, k + 1] -= load * a**2 * (L_elt + 2 * b) / L_elt**3
        self.loads[1, k + 1] += load / L_elt**2 * (b * a**2)

    def add_boundary_moment(self, coordinate: float, value: float):
        """
        adds a moment load onto a node (positive counterclockwise)
        
        :param (float) coordinate: x-coordinate where load is applied. Has to be the first or last node.
        :param (float) load: magnitude of load
        """
        k = np.searchsorted(self.nodes, coordinate)
        if coordinate != self.nodes[k]:
            raise ValueError("Node not in discetization")
        else:
            self.loads[1, k] += value

    def add_boundary_shear(self, coordinate: float, value: float):
        """
        adds a shear load onto a node (positive in y-direction)

        :param (float) coordinate: x-coordinate where load is applied. Has to be the first or last node.
        :param (float) load: magnitude of load
        """
        k = np.searchsorted(self.nodes, coordinate)
        if coordinate != self.nodes[k]:
            raise ValueError("Node not in discetization")
        else:
            self.loads[0, k] += value

    def add_udl(self, interval: list[float], magnitude: float):
        """
        :param (list[float]) interval: x-coordinates of beginning and end point. Have to be in nodes
        :param (float) magnitude: magnitude of the UDL

        adds a distributed load of magnitude on the interval
        """

        xmin, xmax = interval
        if xmin not in self.nodes or xmax not in self.nodes:
            raise ValueError("Interval must be specified between nodes")
        kmin = np.searchsorted(self.nodes, xmin)
        kmax = np.searchsorted(self.nodes, xmax)
        elements = np.arange(kmin, kmax)

        for element in elements:
            L_elt = self.nodes[element + 1] - self.nodes[element]
            self.loads[0, element] -= magnitude * L_elt / 2
            self.loads[0, element + 1] -= magnitude * L_elt / 2
            self.loads[1, element] -= magnitude * L_elt**2 / 12
            self.loads[1, element + 1] += magnitude * L_elt**2 / 12

    def add_triangle_load(self, interval: list[float], magnitude: float):
        """
        Add triangle load on interval

        :param (list[float]) interval: x-coordinates of beginning and end point. Have to be in nodes.
        :param (float) magnitude: magnitude of the load at the beginning of the interval. Load is 0 at its end.
        """
        xmin, xmax = interval
        if xmin not in self.nodes or xmax not in self.nodes:
            raise ValueError("Interval must be specified between nodes")
        kmin = np.searchsorted(self.nodes, xmin)
        kmax = np.searchsorted(self.nodes, xmax)
        if self.nodes[kmin + 1] != self.nodes[kmax]:
            raise ValueError("triangle load can only be applied on one interval")
        L_elt = xmax - xmin
        self.loads[0, kmin] -= 7 * magnitude * L_elt / 20
        self.loads[0, kmax] -= 7 * magnitude * L_elt**2 / 20
        self.loads[1, kmin] -= 3 * magnitude * L_elt / 20
        self.loads[1, kmax] += magnitude * L_elt**2 / 30

    def add_load(self, interval: list[float], w: callable):
        """
        Add load described by function w to the beam. This uses quadrature and thus is much slower that add_udl an. etc.

        :param (list[float]) interval: x-coordinates of beginning and end point. have to be in nodes
        :param (callable) w: integrable function on interval describing load

        adds load described by w on the interval
        """
        xmin, xmax = interval
        if xmin not in self.nodes or xmax not in self.nodes:
            raise ValueError("Interval must be specified between nodes")
        kmin = np.searchsorted(self.nodes, xmin)
        kmax = np.searchsorted(self.nodes, xmax)
        elements = np.arange(kmin, kmax)

        for element in elements:
            xk = self.nodes[element]
            xkp1 = self.nodes[element + 1]

            def func1(x):
                return w(x) * Nk1(x, xk, xkp1)

            def func2(x):
                return w(x) * Nk2(x, xk, xkp1)

            def func3(x):
                return w(x) * Nk3(x, xk, xkp1)

            def func4(x):
                return w(x) * Nk4(x, xk, xkp1)

            self.loads[0, element] += quad(func1, xk, xkp1)[0]
            self.loads[1, element] += quad(func2, xk, xkp1)[0]
            self.loads[0, element + 1] += quad(func3, xk, xkp1)[0]
            self.loads[1, element + 1] += quad(func4, xk, xkp1)[0]

    def stiffness_matrix(self):
        """
        calculate stiffness matrix
        """
        K = np.zeros((2 * self.nnodes, 2 * self.nnodes))
        for k in range(self.nnodes - 1):
            hk = self.nodes[k + 1] - self.nodes[k]
            K[2 * k : 2 * k + 4, 2 * k : 2 * k + 4] += Kloc(hk, self.E, self.I)
        return K

    def fem(self, return_condition_number=False):
        """
        Performs Finite Element method on the beam given the loads and supports and returns the deflection

        :param (bool, optional) return_condition_number: return condition number of stiffness matrix

        :returns (np.ndarray) deflection: Array of deflections
        :returns (float, optional) cond: condition number of stiffness matrix
        """
        free_idx = np.nonzero((self.free).flatten(order="F"))[0]
        zero_idx = np.nonzero((1 - self.free).flatten(order="F"))[0]
        stiffness = self.stiffness_matrix()
        forces = self.loads.flatten(order="F")
        deflection = np.zeros_like(forces)
        K, F = zero_bc(stiffness, forces, zero_idx)
        if return_condition_number:
            cond = np.linalg.cond(K)
        if free_idx.size > 0:
            if np.linalg.matrix_rank(K) < free_idx.size:
                raise ValueError("system not solvable: maybe not enough supports")
        K = lower_banded(K)
        u = la.solveh_banded(K, F, lower=True)
        np.put(deflection, free_idx, u)
        if return_condition_number:
            return deflection, cond
        else:
            return deflection

    def evaluate_fem(self, x: float, sol=None) -> float:
        """
        Evaluates the interpolation of the beam deflection at any x in [0,L]
        it is recommended to calculate the fem solution beforehand and provide it in "sol" if the function is called multiple times

        :param (float) x: x-coordinate
        :param (np.ndarray) sol: solution returned by self.fem
        """
        if sol is None:
            fem = self.fem()
        else:
            fem = sol
        k = self.find_element(x)
        len_k = self.nodes[k + 1] - self.nodes[k]
        f1, m1, f2, m2 = fem[2 * k : 2 * k + 4]
        xi = 2 / len_k * (x - self.nodes[k]) - 1  # change of variables onto [-1,1]
        return (
            f1 * N1(xi)
            + f2 * N3(xi)
            + m1 * len_k / 2 * N2(xi)
            + m2 * len_k / 2 * N4(xi)
        )

    def plot(self):
        """
        Plot beam. The x and y-axis dont have units sine they could from use case to use case.
        """
        length = self.nodes[-1] - self.nodes[0]
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        xs = np.linspace(0, length, 200)
        fem = self.fem()
        ax.plot(xs, [self.evaluate_fem(x, sol=fem) for x in xs], color="black")
        ax.set_ylabel("deflection")
        ax.set_xlabel("position along the beam")
        fig.suptitle(
            "FEM approximation of beam shape with {} nodes".format(self.nnodes)
        )
        ax.set_title(r"E = {}, I = {}".format(self.E, self.I))
        plt.show()