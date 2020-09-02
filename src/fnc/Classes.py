# coding=UTF-8

import numpy as np

class ClosedLoopData():
    """Object collecting closed loop data points
    Attributes:
        updateInitialConditions: function which updates initial conditions and clear the memory
    """
    def __init__(self, dt, Time, v0):
        """Initialization
        Arguments:
            dt: discretization time
            Time: maximum time [s] which can be recorded
            v0: velocity initial condition
        """
        self.dt = dt
        self.Points = int(Time / dt)  # Number of points in the simulation
        self.u = np.zeros((self.Points, 2))  # Initialize the input vector
        self.x = np.zeros((self.Points + 1, 6))  # Initialize state vector (In curvilinear abscissas)
        self.x_glob = np.zeros((self.Points + 1, 6))  # Initialize the state vector in absolute reference frame
        self.SimTime = 0.0
        self.x[0,0] = v0
        self.x_glob[0,0] = v0

    def updateInitialConditions(self, x, x_glob):
        """Clears memory and resets initial condition
        x: initial condition is the curvilinear reference frame
        x_glob: initial condition in the inertial reference frame
        """
        self.x[0, :] = x
        self.x_glob[0, :] = x_glob

        self.x[1:, :] = 0*self.x[1:, :]
        self.x_glob[1:, :] = 0*self.x_glob[1:, :]


class LMPCprediction():
    """Object collecting the predictions and SS at eath time step
    """
    def __init__(self, N, n, d, TimeLMPC, numSS_Points, Laps):
        """
        Initialization:
            N: horizon length
            n, d: input and state dimensions
            TimeLMPC: maximum simulation time length [s] 最大的仿真时间长度
            num_SSpoints: number used to buils SS at each time step 每一步开始时ss集的点个数num_sspoints
        """
        self.PredictedStates = np.zeros((N+1, n, TimeLMPC, Laps))
        self.PredictedInputs = np.zeros((N, d, TimeLMPC, Laps))

        self.SSused   = np.zeros((n , numSS_Points, TimeLMPC, Laps))
        self.Qfunused = np.zeros((numSS_Points, TimeLMPC, Laps))     #在第numss—points个点时，状态（输入）变量里的第timelmpc项，共有laps圈轨迹

class  TLMPCprediction():
    """Object collecting the predictions and SS at eath time step
    """
    def __init__(self, N, n, d, TimeTLMPC, numTSS_Points, TLaps,numTSR_Points):
        """
        Initialization:
            N: horizon length
            n, d: input and state dimensions
            TimeTLMPC: maximum simulation time length [s] 最大的仿真时间长度（状态变量的个数或输入变量的个数）
            num_TSSpoints: number used to buils SS at each time step 每一步开始时ss集的点个数num_sspoints
        """
        self.PredictedStates = np.zeros((N+1, n, TimeTLMPC, TLaps))
        self.PredictedInputs = np.zeros((N, d, TimeTLMPC, TLaps))


        self.Pfunused = np.zeros((numTSS_Points, TimeTLMPC, TLaps))     #在第numss—points个点时，状态（输入）变量里的第timelmpc项，共有laps圈轨迹
        # 相似路径、被采纳的相似轨迹
        self.TSRused = np.zeros((n, numTSR_Points, TimeTLMPC, TLaps))
     #   self.TSRSused = np.zeros((n, numTSS_Points, TimeTLMPC, TLaps))
        self.TSSused = np.zeros((n, numTSS_Points, TimeTLMPC, TLaps))
