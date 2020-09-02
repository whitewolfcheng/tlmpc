#coding=utf-8
import numpy as np
import pdb
import datetime
from Utilities import Curvature, getAngle

class testSimulator():
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
        #for i in range(0, int(ClosedLoopData.Points)):  #在[0,points)中遍历， self.Points = int(Time / dt)，表示仿真中的点个数
        for i in range(0, 1):
            Controller.solve(x[i, :])                   #x中第一个量为预测个数，第二个量为6个状态变量

          #  u[i, :] = Controller.uPred[0,:]
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
          #  x[i + 1, :], x_glob[i + 1, :] = _DynModel(x[i, :], x_glob[i, :], u[i, :], np, ClosedLoopData.dt, self.map.PointAndTangent)
          #  print x[i+1,:]
          #  SimulationTime = i + 1

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
                Controller.addPoint(x[i, :], u[i, :])
            if self.flagTLMPC == 1:
                Controller.addPoint(x[i, :], u[i, :])

            if (self.laps == 1) and (int(np.floor(x[i+1, 4] / (self.map.TrackLength))))>0:
                print "Simulation terminated: Lap completed"
                break

        ClosedLoopData.SimTime = SimulationTime
        print "Number of laps completed: ", int(np.floor(x[-1, 4] / (self.map.TrackLength)))   #floor为向下取整，比如floor（2.1）=2.-1表示取x的最后一次预测