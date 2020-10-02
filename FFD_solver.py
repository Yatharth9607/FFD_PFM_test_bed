# Jacobi solver for Project iterations

import numpy as np
import math
import time

# Define constants

rho = 1.19  # Density of air (kg/m^3)
Cp = 1005  # Specific heat of air (J/kg-K)
visc_k = 1.5462e-5  # Kinematic viscosity of air (m^2/s)
visc_d = 1.84e-5  # Dynamic viscosity of air (N-s/m^2)
k_air = 0.0261  # Thermal conductivity of air (W/m-K)
g = 9.81  # Gravitational acceleration (m/s^2)
beta = 0.00341  # Expansion of coefficient of air (1/K)


class CFD:
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
        self.CFD_parameter = self.solver_config[3]
        self.time_mult = self.CFD_parameter[1]

        # monitor point
        self.monitor_point = self.solver_config[4]
        self.monitor_x = self.monitor_point[0]
        self.monitor_y = self.monitor_point[1]
        self.monitor_U = []
        self.monitor_V = []
        self.monitor_P = []
        self.monitor_T = []

        # Other Variables
        self.U = None
        self.V = None
        self.T = None
        self.P = None

    def initialize(self):

        b = self.opening[0]

        if b == 1:
            f = 0
        else:
            f = (1 / b ** 2) * (1 + 0.5 * (1 - b) ** 0.75 + 1.414 * (1 - b) ** 0.375)

        N = self.solver_config[0]

        dx = 1 / N

        V_1 = self.input_1[0]
        V_2 = self.input_2[0]

        V_left = self.left_wall[0]
        V_right = self.right_wall[0]
        U_top = self.top_wall[0]
        U_bottom = self.bottom_wall[0]

        V_scale = max(
            abs(V_1), abs(V_2), abs(V_left), abs(V_right), abs(U_top), abs(U_bottom), 1
        )

        time = self.CFD_parameter[0]

        dt = self.time_mult / (N * V_scale)

        Nt = int(time / dt)

        a = visc_k * dt / (dx ** 2)

        P_scale = abs(rho * V_scale ** 2) / 2

        T_1 = self.input_1[1]
        T_2 = self.input_2[1]
        T_ref = self.temp[1]

        DT_scale = max(abs(T_1 - T_2), abs(T_1 - T_ref), abs(T_2 - T_ref), 1)

        q = self.q_total / N

        self.U = np.zeros((Nt + 1, N + 1, N + 2))
        self.V = np.zeros((Nt + 1, N + 2, N + 1))
        self.P = np.zeros((N + 2, N + 2))
        self.T = T_ref * np.ones((Nt + 1, N + 2, N + 2))

        self.monitor_U.append(0)
        self.monitor_V.append(0)
        self.monitor_P.append(0)
        self.monitor_T.append(0)

        return f, N, dx, dt, Nt, a, V_scale, P_scale, DT_scale, q

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

        V_out = (V_1 + V_2) / 2

        # left wall velocity and temperature
        self.U[:, 0, :] = 0
        self.T[:, 0, :] = T_left
        self.V[:, 0, :] = V_left

        # right wall velocity and temperature
        self.U[:, N, :] = 0
        self.T[:, N + 1, :] = T_right
        self.V[:, N + 1, :] = V_right

        # top wall velocity and temperature
        self.V[:, :, N] = 0
        self.T[:, :, N + 1] = T_top
        if V_out != 0:
            for i in range(int((3 / 8) * N + 1), int((5 / 8) * N) + 1):
                self.V[:, i, N] = V_out
        self.U[:, :, N + 1] = U_top

        # bottom wall velocity and temperature
        self.V[:, :, 0] = 0
        self.T[:, :, 0] = T_bottom
        if V_1 != 0:
            for i in range(int(N / 8) + 1, int(N / 4) + 1):
                self.V[:, i, 0] = V_1
                self.T[:, i, 0] = T_1
        if V_2 != 0:
            for i in range(int(3 / 4 * N) + 1, int(7 / 8 * N) + 1):
                self.V[:, i, 0] = V_2
                self.T[:, i, 0] = T_2
        self.U[:, :, 0] = U_bottom
        if V_1 != 0:
            for i in range(int(N / 8), int(N / 4) + 1):
                self.U[:, i, 0] = V_1 / (math.tan(theta_1 * 0.017453293))
        if V_2 != 0:
            for i in range(int(3 / 4 * N), int(7 / 8 * N) + 1):
                self.U[:, i, 0] = V_2 / (math.tan(theta_2 * 0.017453293))

    def set_faces(self, N, dx, f):

        FX = np.ones((N + 1, N + 1))
        FY = np.ones((N + 1, N + 1))
        RTX = np.zeros((N + 1, N + 1))
        RTY = np.zeros((N + 1, N + 1))
        RUX = dx / visc_d * np.ones((N, N + 1))
        RUY = dx / visc_d * np.ones((N, N + 1))
        RVX = dx / visc_d * np.ones((N + 1, N))
        RVY = dx / visc_d * np.ones((N + 1, N))
        LossCoeffY = np.zeros((N + 1, N))

        h_left = self.left_wall[2]
        h_right = self.right_wall[2]
        h_top = self.top_wall[2]
        h_bottom = self.bottom_wall[2]

        L = self.opening[1]

        # interior faces and thermal resistance

        RTX[:, :] = dx / k_air
        RTY[:, :] = dx / k_air

        # left and right faces and thermal resistance

        FX[0, :] = 0
        RTX[0, :] = 1 / h_left
        FX[N, :] = 0
        RTX[N, :] = 1 / h_right

        # top and bottom faces and thermal resistance

        FY[:, 0] = 0
        RTY[:, 0] = 1 / h_bottom
        FY[:, N] = 0
        RTY[:, N] = 1 / h_top

        # U diffusion resistance on top and bottom

        RUY[:, 0] = 2 * dx / visc_d
        RUY[:, N] = 2 * dx / visc_d

        # V diffusion resistance on left and right

        RVX[0, :] = 2 * dx / visc_d
        RVX[N, :] = 2 * dx / visc_d

        # Set Loss Coefficient for resistance in Y direction

        LossCoeffY[:, int(N * L)] = f

        return FX, FY, RTX, RTY, RUX, RUY, RVX, RVY, LossCoeffY

    def Solve_CFD(self):

        ref = time.time()

        f, N, dx, dt, Nt, a, V_scale, P_scale, DT_scale, q = self.initialize()

        self.set_bocos(N)

        FX, FY, RTX, RTY, RUX, RUY, RVX, RVY, LossCoeffY = self.set_faces(N, dx, f)

        for t in range(1, Nt + 1):

            if t != 1:
                self.U[t, :, :] = self.U[t - 1, :, :]
                self.U[t - 1, :, :] = self.U[t - 2, :, :]

                self.V[t, :, :] = self.V[t - 1, :, :]
                self.V[t - 1, :, :] = self.V[t - 2, :, :]

                self.T[t, :, :] = self.T[t - 1, :, :]
                self.T[t - 1, :, :] = self.T[t - 2, :, :]

            self.swapU(t)
            self.solveU(N, t, RUX, RUY, dt, dx)

            self.swapV(t)
            self.solveV(N, t, RVX, RVY, LossCoeffY, dt, dx)

            self.swapU(t)
            self.swapV(t)
            self.project(N, t, FX, FY, dt, dx)

            if self.solve_temp == "Y":
                self.swapT(t)
                self.solveT(N, t, q, RTX, RTY, dt, dx)

            self.monitor(V_scale, P_scale, DT_scale, N, t)

        U_col, V_col, T_col, P_col = self.output(N, t)

        timestamp = round(time.time() - ref, 2)

        monitor_data = []
        monitor_data.append(dt * np.arange(Nt + 1))
        monitor_data.append(self.monitor_U)
        monitor_data.append(self.monitor_V)
        monitor_data.append(self.monitor_P)
        monitor_data.append(self.monitor_T)
        # Append one more array to accommodate temperature iterations (T_iter for PFM is different)
        monitor_data.append(dt * np.arange(Nt + 1))

        mass = self.calcMassBalance(N, t)

        return U_col, V_col, T_col, P_col, timestamp, monitor_data, mass

    def swapU(self, t):
        z = self.U[t - 1, :, :]
        self.U[t - 1, :, :] = self.U[t, :, :]
        self.U[t, :, :] = z

    def swapV(self, t):
        z = self.V[t - 1, :, :]
        self.V[t - 1, :, :] = self.V[t, :, :]
        self.V[t, :, :] = z

    def swapT(self, t):
        z = self.T[t - 1, :, :]
        self.T[t - 1, :, :] = self.T[t, :, :]
        self.T[t, :, :] = z

    def solveU(self, N, t, RUX, RUY, dt, dx):
        vel_iter = self.CFD_parameter[2]

        for k in range(1, vel_iter + 1):
            for i in range(1, N):
                for j in range(1, N + 1):
                    UR = (self.U[t - 1, i, j] + self.U[t - 1, i + 1, j]) / 2
                    UL = (self.U[t - 1, i - 1, j] + self.U[t - 1, i, j]) / 2
                    VT = (self.V[t - 1, i, j] + self.V[t - 1, i + 1, j]) / 2
                    VB = (self.V[t - 1, i, j - 1] + self.V[t - 1, i + 1, j - 1]) / 2

                    b = (
                        max(UL, 0) * self.U[t, i - 1, j]
                        + max(VB, 0) * self.U[t, i, j - 1]
                        - min(UR, 0) * self.U[t, i + 1, j]
                        - min(VT, 0) * self.U[t, i, j + 1]
                    )
                    c = max(UL, 0) + max(VB, 0) - min(UR, 0) - min(VT, 0)

                    RUSum = (
                        self.U[t, i, j + 1] / RUY[i, j]
                        + self.U[t, i, j - 1] / RUY[i, j - 1]
                        + self.U[t, i - 1, j] / RUX[i - 1, j]
                        + self.U[t, i + 1, j] / RUX[i, j]
                    )
                    RUSumCoeff = (
                        1 / RUY[i, j]
                        + 1 / RUY[i, j - 1]
                        + 1 / RUX[i - 1, j]
                        + 1 / RUX[i, j]
                    )

                    Num = (
                        self.U[t - 1, i, j] + (dt / dx) * b + (dt / (rho * dx)) * RUSum
                    )
                    Denom = 1 + (dt / dx) * c + (dt / (rho * dx)) * RUSumCoeff

                    self.U[t, i, j] = Num / Denom

    def solveV(self, N, t, RVX, RVY, LossCoeffY, dt, dx):
        vel_iter = self.CFD_parameter[2]
        T_ref = self.temp[1]

        for k in range(1, vel_iter + 1):
            for i in range(1, N + 1):
                for j in range(1, N):
                    UR = (self.U[t - 1, i, j] + self.U[t - 1, i, j + 1]) / 2
                    UL = (self.U[t - 1, i - 1, j] + self.U[t - 1, i - 1, j + 1]) / 2
                    VT = (self.V[t - 1, i, j] + self.V[t - 1, i, j + 1]) / 2
                    VB = (self.V[t - 1, i, j - 1] + self.V[t - 1, i, j]) / 2

                    b = (
                        max(UL, 0) * self.V[t, i - 1, j]
                        + max(VB, 0) * self.V[t, i, j - 1]
                        - min(UR, 0) * self.V[t, i + 1, j]
                        - min(VT, 0) * self.V[t, i, j + 1]
                    )
                    c = max(UL, 0) + max(VB, 0) - min(UR, 0) - min(VT, 0)

                    RVSum = (
                        self.V[t, i, j + 1] / RVY[i, j]
                        + self.V[t, i, j - 1] / RVY[i, j - 1]
                        + self.V[t, i - 1, j] / RVX[i - 1, j]
                        + self.V[t, i + 1, j] / RVX[i, j]
                    )
                    RVSumCoeff = (
                        1 / RVY[i, j]
                        + 1 / RVY[i, j - 1]
                        + 1 / RVX[i - 1, j]
                        + 1 / RVX[i, j]
                    )

                    Force = (
                        -np.sign(self.V[t - 1, i, j])
                        * ((LossCoeffY[i, j] * dt) / (2 * dx))
                        * self.V[t - 1, i, j] ** 2
                        + g * beta * (self.T[t, i, j] - T_ref) * dt
                    )

                    Num = (
                        self.V[t - 1, i, j]
                        + Force
                        + (dt / dx) * b
                        + (dt / (rho * dx)) * RVSum
                    )

                    Denom = 1 + (dt / dx) * c + (dt / (rho * dx)) * RVSumCoeff

                    self.V[t, i, j] = Num / Denom

    def project(self, N, t, FX, FY, dt, dx):
        proj_iter = self.CFD_parameter[4]
        omega = self.solver_config[1]

        div = np.zeros((N + 1, N + 1))

        # Compute divergence at each scalar cell

        for i in range(1, N + 1):
            for j in range(1, N + 1):
                div[i, j] = (1 / dx) * (
                    self.U[t - 1, i, j]
                    - self.U[t - 1, i - 1, j]
                    + self.V[t - 1, i, j]
                    - self.V[t - 1, i, j - 1]
                )

        # Solve for P at each scalar cell using Jacobi solver
        PFSum = np.zeros((N + 1, N + 1))
        for k in range(1, proj_iter + 1):
            P = self.P
            for i in range(1, N + 1):
                for j in range(1, N + 1):

                    PFSum[i, j] = (
                        P[i + 1, j] * FX[i, j]
                        + P[i - 1, j] * FX[i - 1, j]
                        + P[i, j + 1] * FY[i, j]
                        + P[i, j - 1] * FY[i, j - 1]
                    )
                    PFSumCoeff = FX[i, j] + FX[i - 1, j] + FY[i, j] + FY[i, j - 1]
                    self.P[i, j] = (1 - omega) * P[i, j] + omega * (
                        PFSum[i, j] - (rho * dx ** 2 / dt) * div[i, j]
                    ) / PFSumCoeff

        # Compute new U velocities for interior boundaries

        for i in range(1, N):
            for j in range(1, N + 1):
                self.U[t, i, j] = self.U[t - 1, i, j] - (dt / (rho * dx)) * (
                    self.P[i + 1, j] - self.P[i, j]
                )

        # Compute new V velocities for interior boundaries

        for i in range(1, N + 1):
            for j in range(1, N):
                self.V[t, i, j] = self.V[t - 1, i, j] - (dt / (rho * dx)) * (
                    self.P[i, j + 1] - self.P[i, j]
                )

    def solveT(self, N, t, q, RTX, RTY, dt, dx):
        T_iter = self.CFD_parameter[3]

        for k in range(1, T_iter + 1):
            for i in range(1, N + 1):
                for j in range(1, N + 1):
                    UR = self.U[t, i, j]
                    UL = self.U[t, i - 1, j]
                    VT = self.V[t, i, j]
                    VB = self.V[t, i, j - 1]

                    b = (
                        max(UL, 0) * self.T[t, i - 1, j]
                        + max(VB, 0) * self.T[t, i, j - 1]
                        - min(UR, 0) * self.T[t, i + 1, j]
                        - min(VT, 0) * self.T[t, i, j + 1]
                    )
                    c = max(UL, 0) + max(VB, 0) - min(UR, 0) - min(VT, 0)

                    RTSum = (
                        self.T[t, i, j + 1] / RTY[i, j]
                        + self.T[t, i, j - 1] / RTY[i, j - 1]
                        + self.T[t, i - 1, j] / RTX[i - 1, j]
                        + self.T[t, i + 1, j] / RTX[i, j]
                    )
                    RTSumCoeff = (
                        1 / RTY[i, j]
                        + 1 / RTY[i, j - 1]
                        + 1 / RTX[i - 1, j]
                        + 1 / RTX[i, j]
                    )

                    Num = (
                        self.T[t - 1, i, j]
                        + (dt / (rho * Cp * dx ** 2)) * q
                        + (dt / dx) * b
                        + (dt / (rho * Cp * dx)) * RTSum
                    )
                    Denom = 1 + (dt / dx) * c + (dt / (rho * Cp * dx)) * RTSumCoeff

                    self.T[t, i, j] = Num / Denom

    def monitor(self, V_scale, P_scale, DT_scale, N, t):
        T_ref = self.temp[1]
        mon_i = int(N * self.monitor_x)
        mon_j = int(N * self.monitor_y)

        # Save Timestep Data
        self.monitor_U.append(self.U[t, mon_i, mon_j] / V_scale)
        self.monitor_V.append(self.V[t, mon_i, mon_j] / V_scale)
        self.monitor_P.append(self.P[mon_i, mon_j] / P_scale)
        self.monitor_T.append((self.T[t, mon_i, mon_j] - T_ref) / DT_scale)

    def output(self, N, t):
        U_col = np.zeros((N, N))
        V_col = np.zeros((N, N))
        T_col = np.zeros((N, N))
        P_col = np.zeros((N, N))
        # speed = np.zeros((N, N))

        for i in range(1, N + 1):
            for j in range(1, N + 1):
                U_col[i - 1, j - 1] = (self.U[t, i - 1, j] + self.U[t, i, j]) / 2
                V_col[i - 1, j - 1] = (self.V[t, i, j - 1] + self.V[t, i, j]) / 2
                T_col[i - 1, j - 1] = self.T[t, i, j]
                P_col[i - 1, j - 1] = self.P[i, j]

        U_col = np.transpose(U_col)
        V_col = np.transpose(V_col)
        T_col = np.transpose(T_col)
        P_col = np.transpose(P_col)

        # # Compute Speed
        #
        # for i in range(1, N+1):
        #     for j in range(1, N+1):
        #         speed[i-1, j-1] = 0.5 * ((self.U[t,i-1,j] + self.U[t,i,j]) ** 2 + (self.V[t,i,j-1] + self.V[t,i,j]) ** 2) ** 0.5
        return U_col, V_col, T_col, P_col

    def calcMassBalance(self, N, t):
        mass = np.zeros((N, N))
        for i in range(1, N + 1):
            for j in range(1, N + 1):
                mass[i - 1, j - 1] = (
                    self.U[t, i, j]
                    + self.V[t, i, j]
                    - self.U[t, i - 1, j]
                    - self.V[t, i, j - 1]
                )
        return mass
