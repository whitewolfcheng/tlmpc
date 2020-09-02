#coding=utf-8
import numpy as np
import pdb
import datetime
from Utilities import Curvature, getAngle

class Simulator():
    """Vehicle simulator
    Attributes:
        Sim: given a Controller object run closed-loop simulation and write the data in the ClosedLoopData object
    """
    def __init__(self, map, lap = 0, flagLMPC = 0,flagTLMPC=0):
        """Initialization
        map: map
        lap: number of laps to run. If set to 0 then the simulation is completed when ClosedLoopData is full
        flagLMPC: set to 0 for standart controller. Set to 1 for LMPC --> at iteration j add data to SS^{j-1} (look line 9999)
        """
        self.map = map
        self.laps = lap
        self.flagLMPC = flagLMPC
        self.flagTLMPC=flagTLMPC

    def Sim(self, ClosedLoopData, Controller, LMPCprediction=0,TLMPCprediction=0):
        """Simulate closed-loop system
        ClosedLoopData: object where the closed-loop data are written
        Controller: controller used in the closed-loop
        LMPCprediction: object where the open-loop predictions and safe set are stored
        """
        # Assign x = ClosedLoopData.x. IMPORTANT: x and ClosedLoopData.x share the same memory and therefore
        # if you modofy one is like modifying the other one!
        x      = ClosedLoopData.x
        x_glob = ClosedLoopData.x_glob
        u      = ClosedLoopData.u

        SimulationTime = 0
        for i in range(0, int(ClosedLoopData.Points)):  #在[0,points)中遍历， self.Points = int(Time / dt)，表示仿真中的点个数
            Controller.solve(x[i, :])                   #x中第一个量为预测个数，第二个量为6个状态变量

            u[i, :] = Controller.uPred[0,:]
           # print '小车输入'
          #  print u[i,:]

            if LMPCprediction != 0:
                Controller.LapTime = i
                LMPCprediction.PredictedStates[:,:,i, Controller.it]   = Controller.xPred
                LMPCprediction.PredictedInputs[:, :, i, Controller.it] = Controller.uPred
                LMPCprediction.SSused[:, :, i, Controller.it]          = Controller.SS_PointSelectedTot
                LMPCprediction.Qfunused[:, i, Controller.it]           = Controller.Qfun_SelectedTot
            if TLMPCprediction != 0:
                Controller.LapTime = i
                TLMPCprediction.PredictedStates[:,:,i, Controller.it]   = Controller.xPred
                TLMPCprediction.PredictedInputs[:, :, i, Controller.it] = Controller.uPred
                TLMPCprediction.TSRused[:, :, i, Controller.it]         = Controller.TSR_PointSelectedTot
                TLMPCprediction.TSSused[:, :, i, Controller.it]          = Controller.TSS_PointSelectedTot
                TLMPCprediction.Pfunused[:, i, Controller.it]           = Controller.Pfun_SelectedTot
          #  print '小车状态更新'
         #   print x[i, :]
            x[i + 1, :], x_glob[i + 1, :] = _DynModel(x[i, :], x_glob[i, :], u[i, :], np, ClosedLoopData.dt, self.map.PointAndTangent)
          #  print x[i+1,:]
            SimulationTime = i + 1

            # print(getAngle(x[i+1,4], x[i+1,3], self.map.PointAndTangent))
            # print(x[i+1,3], x_glob[i+1,3], wrap(x_glob[i+1,3]))
            # pdb.set_trace()

#            if i <= 5:
#                print("Linearization time: %.4fs Solver time: %.4fs" % (Controller.linearizationTime.total_seconds(), Controller.solverTime.total_seconds()))
#                print "Time: ", i * ClosedLoopData.dt, "Current State and Input: ", x[i, :], u[i, :]

            if Controller.feasible == 0:
                print "Unfeasible at time ", i*ClosedLoopData.dt
                print "Cur State: ", x[i, :], "Iteration ", Controller.it
                break

            if self.flagLMPC == 1:
                Controller.addPoint(x[i, :], u[i, :],x_glob[i,:])
            if self.flagTLMPC == 1:
                Controller.addPoint(x[i, :], u[i, :],x_glob[i,:])

            if (self.laps == 1) and (int(np.floor(x[i+1, 4] / (self.map.TrackLength))))>0:
                print "Simulation terminated: Lap completed"
                break

        ClosedLoopData.SimTime = SimulationTime
        print "Number of laps completed: ", int(np.floor(x[-1, 4] / (self.map.TrackLength)))   #floor为向下取整，比如floor（2.1）=2.-1表示取x的最后一次预测

class PID:
    """Create the PID controller used for path following at constant speed
    Attributes:
        solve: given x0 computes the control action
    """
    def __init__(self, vt):
        """Initialization
        Arguments:
            vt: target velocity
        """
        self.vt = vt
        self.uPred = np.zeros([1,2])      #输入初始化

        startTimer = datetime.datetime.now()
        endTimer = datetime.datetime.now(); deltaTimer = endTimer - startTimer #datetime函数用来输出当前时刻（年月日时分秒），deltatimer求出间隔时刻
        self.solverTime = deltaTimer
        self.linearizationTime = deltaTimer
        self.feasible = 1

    def solve(self, x0):
        """Computes control action
        Arguments:
            x0: current state position 当前状态
        """
        vt = self.vt
        self.uPred[0, 0] = - 0.6 * x0[5] - 0.9 * x0[3] + np.max([-0.9, np.min([np.random.randn() * 0.25, 0.9])])  #输出角度 其中x[5]为ey，x[3]为epsi， P调节后再加上-0.9到0.9之间的随机误差
        self.uPred[0, 1] = 1.5 * (vt - x0[0]) + np.max([-0.2, np.min([np.random.randn() * 0.10, 0.2])])  #输出加速度 其中x[0]为纵向速度 ，P调节后再加上-0.2到0.2之间的随机误差
# ======================================================================================================================
# ======================================================================================================================
# ================================ Internal functions for change of coordinates ========================================
# ======================================================================================================================
# ======================================================================================================================
def _DynModel(x, x_glob, u, np, dt, PointAndTangent):
    # This function computes the system evolution. Note that the discretization is deltaT and therefore is needed that
    # dt <= deltaT and ( dt / deltaT) = integer value

    # Vehicle Parameters
    m  = 1.98
    lf = 0.125
    lr = 0.125
    Iz = 0.024
    Df = 0.8 * m * 9.81 / 2.0
    Cf = 1.25
    Bf = 1.0
    Dr = 0.8 * m * 9.81 / 2.0
    Cr = 1.25
    Br = 1.0

    # Discretization Parameters
    deltaT = 0.001
    x_next     = np.zeros(x.shape[0])
    cur_x_next = np.zeros(x.shape[0])

    # Extract the value of the states
    delta = u[0]
    a     = u[1]

    psi = x_glob[3]
    X = x_glob[4]
    Y = x_glob[5]

    vx    = x[0]
    vy    = x[1]
    wz    = x[2]
    epsi  = x[3]
    s     = x[4]
    ey    = x[5]

    # Initialize counter 计数器
    i = 0
    while( (i+1) * deltaT <= dt):
        # Compute tire split angle
        alpha_f = delta - np.arctan2( vy + lf * wz, vx )
        alpha_r = - np.arctan2( vy - lf * wz , vx)

        # Compute lateral force at front and rear tire
        Fyf = 2 * Df * np.sin( Cf * np.arctan(Bf * alpha_f ) )
        Fyr = 2 * Dr * np.sin( Cr * np.arctan(Br * alpha_r ) )

        # Propagate the dynamics of deltaT    #直角坐标系
        x_next[0] = vx  + deltaT * (a - 1 / m * Fyf * np.sin(delta) + wz*vy)
        x_next[1] = vy  + deltaT * (1 / m * (Fyf * np.cos(delta) + Fyr) - wz * vx)
        x_next[2] = wz  + deltaT * (1 / Iz *(lf * Fyf * np.cos(delta) - lr * Fyr) )
        x_next[3] = psi + deltaT * (wz)
        x_next[4] =   X + deltaT * ((vx * np.cos(psi) - vy * np.sin(psi)))
        x_next[5] =   Y + deltaT * (vx * np.sin(psi)  + vy * np.cos(psi))

        cur = Curvature(s, PointAndTangent)    #曲率坐标系下
        cur_x_next[0] = vx   + deltaT * (a - 1 / m * Fyf * np.sin(delta) + wz*vy)
        cur_x_next[1] = vy   + deltaT * (1 / m * (Fyf * np.cos(delta) + Fyr) - wz * vx)
        cur_x_next[2] = wz   + deltaT * (1 / Iz *(lf * Fyf * np.cos(delta) - lr * Fyr) )
        cur_x_next[3] = epsi + deltaT * ( wz - (vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey) * cur )
        cur_x_next[4] = s    + deltaT * ( (vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey) )
        cur_x_next[5] = ey   + deltaT * (vx * np.sin(epsi) + vy * np.cos(epsi))

        # Update the value of the states
        psi  = x_next[3]
        X    = x_next[4]
        Y    = x_next[5]

        vx   = cur_x_next[0]
        vy   = cur_x_next[1]
        wz   = cur_x_next[2]
        epsi = cur_x_next[3]
        s    = cur_x_next[4]
        ey   = cur_x_next[5]

        if (s < 0):
            print "Start Point: ", x, " Input: ", u
            print "x_next: ", x_next

        # Increment counter
        i = i+1

    # Noises 应该是人为加入随机噪声，使仿真结果更真实
   # np.random.seed(0)
    noise_vx = np.max([-0.05, np.min([np.random.randn() * 0.01, 0.05])]) #噪声在-0.05到0.05之间随机产生
    noise_vy = np.max([-0.1, np.min([np.random.randn() * 0.01, 0.1])])
    noise_wz = np.max([-0.05, np.min([np.random.randn() * 0.005, 0.05])])

    cur_x_next[0] = cur_x_next[0] + 0.1*noise_vx
    cur_x_next[1] = cur_x_next[1] + 0.1*noise_vy
    cur_x_next[2] = cur_x_next[2] + 0.1*noise_wz

    return cur_x_next, x_next  #曲率坐标系和直角坐标系


def wrap(angle):  #角度大于π就让角度减去2π，小于π就让角度加上2π，保证角度在-π，π之间
    if angle < -np.pi:
        w_angle = 2 * np.pi + angle
    elif angle > np.pi:
        w_angle = angle - 2 * np.pi
    else:
        w_angle = angle

    return w_angle