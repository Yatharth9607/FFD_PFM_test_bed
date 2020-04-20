# Gauss-seidel solver is used

import numpy as np
import math
import time

# Constants

rho = 1.19  # Density of air (kg/m^3)
Cp = 1005  # Specific heat of air (J/kg-K)
visc_k = 1.5462e-5  # Kinematic viscosity of air (m^2/s)
visc_d = 1.84e-5  # Dynamic viscosity of air (N-s/m^2)
k_air = 0.0261  # Thermal conductivity of air (W/m-K)
g = 9.81  # Gravitational acceleration (m/s^2)
beta = 0.00341  # Expansion of coefficient of air (1/K)

class PFM:

    def __init__(self, design_inputs, solver_config):

        # Input - 1 design variables
        self.input_1 = design_inputs[0]

        # Input - 2 design variables
        self.input_2 = design_inputs[1]

        # Heat source
        self.q_total = design_inputs[2]

        # Opening parameters
        self.opening = design_inputs[3]

        # Left wall parameters
        self.left_wall = design_inputs[4]

        # Right wall parameters
        self.right_wall = design_inputs[5]

        # Top wall parameters
        self.top_wall = design_inputs[6]

        # Bottom wall parameters
        self.bottom_wall = design_inputs[7]

        # Solver config
        self.solver_config = solver_config

        # Temperature parameters
        self.temp = self.solver_config[2]
        self.solve_temp = self.temp[0]

        # Iteration parameters
        self.PFM_parameter = self.solver_config[5]

        # monitor point
        self.monitor = self.solver_config[4]
        self.monitor_x = self.monitor[0]
        self.monitor_y = self.monitor[1]
        self.monitor_U = []
        self.monitor_V = []
        self.monitor_P = []
        self.monitor_T = []

        # Other Variables
        self.U = None
        self.V = None
        self.T = None
        self.P = None

        self.U_col = None
        self.V_col = None
        self.T_col = None
        self.P_col = None

    def initialize(self):

        b = self.opening[0]

        if b == 1:
            f = 0
        else:
            f = (1 / b ** 2) * (1 + 0.5 * (1 - b) ** 0.75 + 1.414 * (1 - b) ** 0.375)

        N = self.solver_config[0]

        dx = 1/N

        V_1 = self.input_1[0]
        V_2 = self.input_2[0]

        V_left = self.left_wall[0]
        V_right = self.right_wall[0]
        U_top = self.top_wall[0]
        U_bottom = self.bottom_wall[0]

        V_scale = max(abs(V_1), abs(V_2), abs(V_left), abs(V_right), abs(U_top), abs(U_bottom), 1)

        P_scale = abs(rho*V_scale**2)/2

        T_1 = self.input_1[1]
        T_2 = self.input_2[1]
        T_ref = self.temp[1]

        DT_scale = max(abs(T_1 - T_2), abs(T_1 - T_ref), abs(T_2 - T_ref), 1)

        q = self.q_total/N

        self.U = np.zeros((N+1, N+2))
        self.V = np.zeros((N+2, N+1))
        self.P = np.zeros((N+2, N+2))
        self.T = T_ref*np.ones((N+2, N+2))
        self.phi = np.zeros((N + 2, N + 2))

        self.U_col = np.zeros((N, N))
        self.V_col = np.zeros((N, N))
        self.T_col = np.zeros((N, N))
        self.P_col = np.zeros((N, N))

        self.monitor_U.append(0)
        self.monitor_V.append(0)
        self.monitor_P.append(0)
        self.monitor_T.append(0)

        return f, N, dx, V_scale, P_scale, DT_scale, q

    def set_bocos(self, N):

        # N = self.solver_config[0]

        V_1 = self.input_1[0]
        theta_1 = self.input_1[2]
        T_1 = self.input_1[1]

        V_2 = self.input_2[0]
        theta_2 = self.input_2[2]
        T_2 = self.input_2[1]

        V_left = self.left_wall[0]
        T_left = self.left_wall[1]

        V_right = self.right_wall[0]
        T_right = self.right_wall[1]

        U_top = self.top_wall[0]
        T_top = self.top_wall[1]

        U_bottom = self.bottom_wall[0]
        T_bottom = self.bottom_wall[1]

        V_out = (V_1 + V_2)/2

        # left wall velocity and temperature
        self.U[0,:] = 0
        self.T[0,:] = T_left
        self.V[0,:] = V_left

        # right wall velocity and temperature
        self.U[N,:] = 0
        self.T[N+1,:] = T_right
        self.V[N+1,:] = V_right

        # top wall velocity and temperature
        self.V[:,N] = 0
        self.T[:,N+1] = T_top
        if V_out != 0:
            for i in range(int((3/8)*N+1), int((5/8)*N)+1):
                self.V[i,N] = V_out
        self.U[:,N+1] = U_top

        # bottom wall velocity and temperature
        self.V[:,0] = 0
        self.T[:,0] = T_bottom
        if V_1 != 0:
            for i in range(int(N/8)+1, int(N/4)+1):
                self.V[i,0] = V_1
                self.T[i,0] = T_1
        if V_2 != 0:
            for i in range(int(3/4*N)+1, int(7/8*N)+1):
                self.V[i,0] = V_2
                self.T[i,0] = T_2
        self.U[:,0] = U_bottom
        if V_1 != 0:
            for i in range(int(N/8), int(N/4)+1):
                self.U[i,0] = V_1/(math.tan(theta_1*0.017453293))
        if V_2 != 0:
            for i in range(int(3/4*N), int(7/8*N)+1):
                self.U[i,0] = V_2/(math.tan(theta_2*0.017453293))

    def checkMassBalanceBottom(self, d):
        f = self.opening[0]
        L = self.opening[1]
        N = self.solver_config[0]
        V_1 = self.input_1[0]
        V_2 = self.input_2[0]
        sum = 0
        for i in range(1, N+1):
            sum = sum + np.sign(self.P[i, int(N * L)] - self.P[i, int(N * L) + 1] + d) * (
                    (2 * abs(self.P[i, int(N * L)] - self.P[i, int(N * L) + 1] + d)) / (f * rho)) ** 0.5

        massBalanceBottom = (N / 8) * (V_1 + V_2) - sum

        return massBalanceBottom

    def presCorrectionBottom(self):
        low = -5000
        high = 5000
        lowVal = self.checkMassBalanceBottom(low)
        highVal = self.checkMassBalanceBottom(high)

        for ii in range(1, 50):
            Middle = (low + high) / 2
            midVal = self.checkMassBalanceBottom(Middle)

            if np.sign(midVal) == np.sign(lowVal):
                low = Middle
                lowVal = midVal
            else:
                high = Middle
                highVal = midVal

        pressureCorrectionBottom = Middle

        return pressureCorrectionBottom

    def Solve_PFM(self):

        ref = time.time()

        f, N, dx, V_scale, P_scale, DT_scale, q = self.initialize()

        self.set_bocos(N)

        PFMOuter_Iter = self.PFM_parameter[0]
        T_iter = self.PFM_parameter[2]
        for t in range(1, PFMOuter_Iter+1):
            if f == 0:

                self.velPotAll(N, dx)

                self.BernoulliPres(N)

            else:

                self.BernoulliPres(N)

                self.scalePresBottom(N)

                self.resetMidV(N)

                self.velPotBottom(N, dx)

                self.velPotTop(N, dx)

                self.monitor_UVP(V_scale, P_scale, N)

        if self.solve_temp == 'Y':

            self.temperature(DT_scale, N)

        self.output(N)

        timestamp = round(time.time() - ref, 2)

        monitor_data = []
        monitor_data.append(np.arange(PFMOuter_Iter + 1))
        monitor_data.append(self.monitor_U)
        monitor_data.append(self.monitor_V)
        monitor_data.append(self.monitor_P)
        monitor_data.append(self.monitor_T)
        # Append one more array to accommodate temperature iterations (T_iter for PFM is different)
        monitor_data.append(np.arange(T_iter + 1))

        mass = self.calcMassBalance(N)

        return self.U_col, self.V_col, self.T_col, self.P_col, timestamp, monitor_data, mass

    def velPotAll(self, N, dx):
        # Velocity Potentials All() ---------------------

        VelPot_Iter = self.PFM_parameter[1]
        omega = self.solver_config[1]
        # Solve for velocity potentials at each scalar cell using Gauss-Seidel

        for k in range(1, VelPot_Iter + 1):

            # Bottom and top fictitious halo cells
            for i in range(1, N + 1):
                self.phi[i, 0] = (1 - omega) * self.phi[i, 0] + omega * (self.phi[i, 1] - dx * self.V[i, 0])
                self.phi[i, N + 1] = (1 - omega) * self.phi[i, N + 1] + omega * (self.phi[i, N] + dx * self.V[i, N])

            # left and right fictitious halo cells
            for j in range(1, N + 1):
                self.phi[0, j] = (1 - omega) * self.phi[0, j] + omega * (self.phi[1, j] - dx * self.U[0, j])
                self.phi[N + 1, j] = (1 - omega) * self.phi[N + 1, j] + omega * (self.phi[N, j] + dx * self.U[N, j])

            # interior cells
            for i in range(1, N + 1):
                for j in range(1, N + 1):
                    self.phi[i, j] = (1 - omega) * self.phi[i, j] + (omega / 4) * (self.phi[i + 1, j] + self.phi[i, j + 1] + self.phi[i - 1, j] + self.phi[i, j - 1])

            # compute interior velocities from velocity potentials
            for i in range(1, N):
                for j in range(1, N + 1):
                    self.U[i, j] = (1 / dx) * (self.phi[i + 1, j] - self.phi[i, j])

            for i in range(1, N + 1):
                for j in range(1, N):
                    self.V[i, j] = (1 / dx) * (self.phi[i, j + 1] - self.phi[i, j])

    def BernoulliPres(self, N):
        # BernoulliPressure() --------------------------

        # Interpolate Vectors() -----------------------

        for i in range(1, N + 1):
            for j in range(1, N + 1):
                self.U_col[i - 1, j - 1] = (self.U[i - 1, j] + self.U[i, j]) / 2
                self.V_col[i - 1, j - 1] = (self.V[i, j - 1] + self.V[i, j]) / 2
                self.T_col[i - 1, j - 1] = self.T[i, j]
                self.P_col[i - 1, j - 1] = self.P[i, j]

        # compute relative pressure

        for j in range(N, 0, -1):
            for i in range(N, 0, -1):
                self.P[i, j] = -(rho / 2) * (self.U_col[i-1, j-1] ** 2 + self.V_col[i-1, j-1] ** 2)

    def scalePresBottom(self, N):
        # ScalePressuresBottom() ------------------------

        d = self.presCorrectionBottom()
        L = self.opening[1]
        # Scale pressures with new pressure correction value

        for i in range(1, N + 1):
            for j in range(1, int(N * L) + 1):
                self.P[i, j] = self.P[i, j] + d

    def resetMidV(self, N):
        # ResetMidV() --------------------------------
        # Resets the top velocity BC
        f = self.opening[0]
        L = self.opening[1]
        for i in range(1, N + 1):
            self.V[i, int(N * L)] = np.sign(self.P[i, int(N * L)] - self.P[i, int(N * L) + 1]) * (
                    (2 * abs(self.P[i, int(N * L)] - self.P[i, int(N * L) + 1])) / (f * rho)) ** 0.5

    def velPotBottom(self, N, dx):
        # Velocity_Potentials_Bottom() ---------------------------------

        # Solve for velocity potentials at each scalar cell using Gauss-Seidel
        VelPot_Iter = self.PFM_parameter[1]
        omega = self.solver_config[1]
        L = self.opening[1]

        for k in range(1, VelPot_Iter + 1):     # number of velocity potential iterations

            # bottom and top fictitious halo cells
            for i in range(1, N + 1):
                self.phi[i, 0] = (1 - omega) * self.phi[i, 0] + omega * (self.phi[i, 1] - dx * self.V[i, 0])
                self.phi[i, int(N * L) + 1] = (1 - omega) * self.phi[i, int(N * L) + 1] + omega * (
                            self.phi[i, int(N * L)] + dx * self.V[i, int(N * L)])

            # left and right fictitious halo cells
            for j in range(1, int(N * L) + 1):
                self.phi[0, j] = (1 - omega) * self.phi[0, j] + omega * (self.phi[1, j] - dx * self.U[0, j])
                self.phi[N + 1, j] = (1 - omega) * self.phi[N + 1, j] + omega * (self.phi[N, j] + dx * self.U[N, j])

            # interior cells
            for i in range(1, N + 1):
                for j in range(1, int(N * L) + 1):
                    self.phi[i, j] = (1 - omega) * self.phi[i, j] + (omega / 4) * (
                            self.phi[i + 1, j] + self.phi[i, j + 1] + self.phi[i - 1, j] + self.phi[i, j - 1])

        # compute interior velocities from velocity potentials
        for i in range(1, N):
            for j in range(1, int(N * L) + 1):
                self.U[i, j] = (1 / dx) * (self.phi[i + 1, j] - self.phi[i, j])

        for i in range(1, N + 1):
            for j in range(1, int(N * L)):
                self.V[i, j] = (1 / dx) * (self.phi[i, j + 1] - self.phi[i, j])

    def velPotTop(self, N, dx):
        # Velocity_Potentials_Top() -------------------------------------

        # Solve for velocity potentials at each scalar cell using Gauss-Seidel
        VelPot_Iter = self.PFM_parameter[1]
        L = self.opening[1]
        omega = self.solver_config[1]

        for k in range(1, VelPot_Iter):  # number of velocity potential iterations

            # bottom and top fictitious halo cells
            for i in range(1, N + 1):
                self.phi[i, int(N * L)] = (1 - omega) * self.phi[i, int(N * L)] + omega * (
                            self.phi[i, int(N * L) + 1] - dx * self.V[i, int(N * L)])
                self.phi[i, N + 1] = (1 - omega) * self.phi[i, N + 1] + omega * (self.phi[i, N] + dx * self.V[i, N])

            # left and right fictitious halo cells
            for j in range(int(N * L) + 1, N + 1):
                self.phi[0, j] = (1 - omega) * self.phi[0, j] + omega * (self.phi[1, j] - dx * self.U[0, j])
                self.phi[N + 1, j] = (1 - omega) * self.phi[N + 1, j] + omega * (self.phi[N, j] + dx * self.U[N, j])

            # interior cells
            for i in range(1, N + 1):
                for j in range(int(N * L) + 1, N + 1):
                    self.phi[i, j] = (1 - omega) * self.phi[i, j] + (omega / 4) * (
                                self.phi[i + 1, j] + self.phi[i, j + 1] + self.phi[i - 1, j] + self.phi[i, j - 1])

        # compute interior velocities from velocity potentials

        for i in range(1, N):
            for j in range(int(N * L) + 1, N+1):
                self.U[i, j] = (1 / dx) * (self.phi[i + 1, j] - self.phi[i, j])

        for i in range(1, N + 1):
            for j in range(int(N * L) + 1, N):
                self.V[i, j] = (1 / dx) * (self.phi[i, j + 1] - self.phi[i, j])

    def temperature(self, DT_scale, N):
        # Temperature() -------------------------------
        Temp_Iter = self.PFM_parameter[2]

        for k in range(1, Temp_Iter + 1):
            for i in range(1, N + 1):
                for j in range(1, N + 1):
                    b = max(self.U[i - 1, j], 0) * self.T[i - 1, j] + max(self.V[i, j - 1], 0) * self.T[i, j - 1] - min(self.U[i, j], 0) * self.T[
                        i + 1, j] - min(self.V[i, j], 0) * self.T[i, j + 1]
                    c = max(self.U[i, j], 0) + max(self.V[i, j], 0) - min(self.U[i - 1, j], 0) - min(self.V[i, j - 1], 0)
                    self.T[i, j] = b / c

            self.monitor_temp(DT_scale, N)

    def monitor_UVP(self, V_scale, P_scale, N):
        mon_i = int(N * self.monitor_x)
        mon_j = int(N * self.monitor_y)

        # Save Timestep Data
        self.monitor_U.append(self.U[mon_i, mon_j] / V_scale)
        self.monitor_V.append(self.V[mon_i, mon_j] / V_scale)
        self.monitor_P.append(self.P[mon_i, mon_j] / P_scale)

    def monitor_temp(self, DT_scale, N):
        T_ref = self.temp[1]
        mon_i = int(N * self.monitor_x)
        mon_j = int(N * self.monitor_y)

        # Save timestep data
        self.monitor_T.append((self.T[mon_i, mon_j] - T_ref)/ DT_scale)

    def output(self, N):
        for i in range(1, N+1):
            for j in range(1, N+1):
                self.U_col[i-1, j-1] = (self.U[i-1,j] + self.U[i,j]) / 2
                self.V_col[i-1, j-1] = (self.V[i,j-1] + self.V[i,j]) / 2
                self.T_col[i-1, j-1] = self.T[i,j]
                self.P_col[i-1, j-1] = self.P[i,j]

        self.U_col = np.transpose(self.U_col)
        self.V_col = np.transpose(self.V_col)
        self.T_col = np.transpose(self.T_col)
        self.P_col = np.transpose(self.P_col)

    def calcMassBalance(self, N):
        mass = np.zeros((N, N))
        for i in range(1, N + 1):
            for j in range(1, N + 1):
                mass[i - 1, j - 1] = self.U[i, j] + self.V[i, j] - self.U[i - 1, j] - self.V[i, j - 1]
        return mass