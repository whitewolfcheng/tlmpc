
# coding=UTF-8
# ----------------------------------------------------------------------------------------------------------------------
# Licensing Information: You are free to use or extend these projects for
# education or reserach purposes provided that you provide clear attribution to UC Berkeley,
# including a reference to the papers describing the control framework:
# [1] Ugo Rosolia and Francesco Borrelli. "Learning Model Predictive Control for Iterative Tasks. A Data-Driven
#     Control Framework." In IEEE Transactions on Automatic Control (2017).
#
# [2] Ugo Rosolia, Ashwin Carvalho, and Francesco Borrelli. "Autonomous racing using learning model predictive control."
#     In 2017 IEEE American Control Conference (ACC)
#
# [3] Maximilian Brunner, Ugo Rosolia, Jon Gonzales and Francesco Borrelli "Repetitive learning model predictive
#     control: An autonomous racing example" In 2017 IEEE Conference on Decision and Control (CDC)
#
# [4] Ugo Rosolia and Francesco Borrelli. "Learning Model Predictive Control for Iterative Tasks: A Computationally
#     Efficient Approach for Linear System." IFAC-PapersOnLine 50.1 (2017).
#
# Attibution Information: Code developed by Ugo Rosolia
# (for clarifiactions and suggestions please write to ugo.rosolia@berkeley.edu).
#
# Code description: Simulation of the Learning Model Predictive Controller (LMPC). The main file runs:
# 1) A PID path following controller
# 2) A Model Predictive Controller (MPC) which uses a LTI model identified from the data collected with the PID in 1)
# 3) A MPC which uses a LTV model identified from the date collected in 1)
# 4) A LMPC for racing where the safe set and value function approximation are build using the data from 1), 2) and 3)
# ----------------------------------------------------------------------------------------------------------------------
import sys
sys.path.append('fnc')
from SysModel import Simulator, PID
from Classes import ClosedLoopData, LMPCprediction, TLMPCprediction,ClosedLoopDataLMPC
from PathFollowingLTVMPC import PathFollowingLTV_MPC
from PathFollowingLTIMPC import PathFollowingLTI_MPC
from Track import Map, unityTestChangeOfCoordinates
from LMPC import ControllerLMPC
from TransferLMPC import ControllerTLMPC
from Utilities import Regression
from plot import plotTrajectory, plotClosedLoopLMPC, animation_xy, animation_states, saveGif_xyResults, Save_statesAnimation,TLMPC_xyResults
import numpy as np
import matplotlib.pyplot as plt
import pdb
import  pickle
from testsysmodel import testSimulator

# ======================================================================================================================
# ============================ Choose which controller to run ==========================================================
# ======================================================================================================================
RunPID     = 0; plotFlag       = 0
RunMPC     = 0; plotFlagMPC    = 0
RunMPC_tv  = 0; plotFlagMPC_tv = 0
RunLMPC    = 0; plotFlagLMPC   = 0; animation_xyFlag = 0; animation_stateFlag = 0
RunTLMPC    = 1; plotFlagTLMPC   = 1
# ======================================================================================================================
# ============================ Initialize parameters for path following ================================================
# ======================================================================================================================
dt         = 1.0/10.0        # Controller discretization time
Time       = 100             # Simulation time for PID
TimeMPC    = 100             # Simulation time for path following MPC
TimeMPC_tv = 100             # Simulation time for path following LTV-MPC
vt         = 0.8             # Reference velocity for path following controllers
v0         = 0.5             # Initial velocity at lap 0
N          = 12              # Horizon length
n = 6;   d = 2               # State and Input dimension

# Path Following tuning
Q = np.diag([1.0, 1.0, 1, 1, 0.0, 100.0]) # vx, vy, wz, epsi, s, ey
R = np.diag([1.0, 10.0])                  # delta, a

map = Map(0.4)                            # Initialize the map
simulator = Simulator(map)                # Initialize the Simulator

# ======================================================================================================================
# ==================================== Initialize parameters for LMPC ==================================================
# ======================================================================================================================
TimeLMPC   = 400              # Simulation time
Laps       = 40+2           # Total LMPC laps

# Safe Set Parameters
LMPC_Solver = "CVX"           # Can pick CVX for cvxopt or OSQP. For OSQP uncomment line 14 in LMPC.py
numSS_it = 2                  # Number of trajectories used at each iteration to build the safe set
numSS_Points = 40         # Number of points to select from each trajectory to build the safe set

# Tuning Parameters
Qslack  =  2 * 5 * np.diag([10, 1, 1, 1, 10, 1])            # Cost on the slack variable for the terminal constraint
Qlane   =  1 * np.array([0, 10])                        # Quadratic and linear slack lane cost
Q_LMPC  =  0 * np.diag([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # State cost x = [vx, vy, wz, epsi, s, ey]
R_LMPC  =  0 * np.diag([1.0, 1.0])                      # Input cost u = [delta, a]
dR_LMPC = 10 * np.array([1.0, 10.0])                    # Input rate cost u
#输入约束的边界值
inputConstr = np.array([[0.5, 0.5],                     # Min Steering and Max Steering
                        [1.0, 1.0]])                    # Min Acceleration and Max Acceleration

# Initialize LMPC simulator
LMPCSimulator = Simulator(map, 1, 1)

# ======================================================================================================================
# ==================================== Initialize parameters for TLMPC ==================================================
# ======================================================================================================================
TimeTLMPC   = 400   # Simulation time
TLaps       = 1           # Total LMPC laps
# Safe Set Parameters
TLMPC_Solver = "CVX"           # Can pick CVX for cvxopt or OSQP. For OSQP uncomment line 14 in LMPC.py
numTSS_it = 2                  # Number of trajectories used at each iteration to build the safe set
numTSS_Points = 40         # Number of points to select from each trajectory to build the safe set
numSR_Points=30   #迁移学习中相似路径的片段点个数
# Tuning Parameters
TQslack  =  2 * 5 * np.diag([10, 1, 1, 1, 10, 1])            # Cost on the slack variable for the terminal constraint
TQlane   =  1 * np.array([0, 10])                        # Quadratic and linear slack lane cost
Q_TLMPC  =  0 * np.diag([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # State cost x = [vx, vy, wz, epsi, s, ey]
R_TLMPC  =  0 * np.diag([1.0, 1.0])                      # Input cost u = [delta, a]
dR_TLMPC = 10 * np.array([1.0, 10.0])                    # Input rate cost u
#输入约束的边界值
inputConstr = np.array([[0.5, 0.5],                     # Min Steering and Max Steering
                        [1.0, 1.0]])                    # Min Acceleration and Max Acceleration

# Initialize LMPC simulator
TLMPCSimulator = Simulator(map, 1, 0,1)
# ======================================================================================================================
# ======================================= PID path following ===========================================================
# ======================================================================================================================
print "Starting PID"
if RunPID == 1:
    ClosedLoopDataPID = ClosedLoopData(dt, Time , v0)
    PIDController = PID(vt)
    simulator.Sim(ClosedLoopDataPID, PIDController)

    file_data = open(sys.path[0]+'/data/ClosedLoopDataPID.obj', 'wb')
    pickle.dump(ClosedLoopDataPID, file_data)
    file_data.close()
else:
    file_data = open(sys.path[0]+'/data/ClosedLoopDataPID.obj', 'rb')
    ClosedLoopDataPID = pickle.load(file_data)
    file_data.close()
print "===== PID terminated"

# ======================================================================================================================
# ======================================  LINEAR REGRESSION ============================================================
# ======================================================================================================================
print "Starting MPC"
lamb = 0.0000001 #应该是增量
A, B, Error = Regression(ClosedLoopDataPID.x, ClosedLoopDataPID.u, lamb)

if RunMPC == 1:
    ClosedLoopDataLTI_MPC = ClosedLoopData(dt, TimeMPC, v0)
    Controller_PathFollowingLTI_MPC = PathFollowingLTI_MPC(A, B, Q, R, N, vt, inputConstr)
    simulator.Sim(ClosedLoopDataLTI_MPC, Controller_PathFollowingLTI_MPC)

    file_data = open(sys.path[0]+'/data/ClosedLoopDataLTI_MPC.obj', 'w')
    pickle.dump(ClosedLoopDataLTI_MPC, file_data)
    file_data.close()
else:
    file_data = open(sys.path[0]+'/data/ClosedLoopDataLTI_MPC.obj', 'r')
    ClosedLoopDataLTI_MPC = pickle.load(file_data)
    file_data.close()
print "===== MPC terminated"

# ======================================================================================================================
# ===================================  LOCAL LINEAR REGRESSION =========================================================
# ======================================================================================================================
print "Starting TV-MPC"
if RunMPC_tv == 1:
    ClosedLoopDataLTV_MPC = ClosedLoopData(dt, TimeMPC_tv, v0)
    Controller_PathFollowingLTV_MPC = PathFollowingLTV_MPC(Q, R, N, vt, n, d, ClosedLoopDataPID.x, ClosedLoopDataPID.u, dt, map, inputConstr)
    simulator.Sim(ClosedLoopDataLTV_MPC, Controller_PathFollowingLTV_MPC)

    file_data = open(sys.path[0]+'/data/ClosedLoopDataLTV_MPC.obj', 'w')
    pickle.dump(ClosedLoopDataLTV_MPC, file_data)
    file_data.close()
else:
    file_data = open(sys.path[0]+'/data/ClosedLoopDataLTV_MPC.obj', 'r')
    ClosedLoopDataLTV_MPC = pickle.load(file_data)
    file_data.close()
print "===== TV-MPC terminated"
# ======================================================================================================================
# ==============================  LMPC w\ LOCAL LINEAR REGRESSION ======================================================
# ======================================================================================================================
print "Starting LMPC"
ClosedLoopLMPC = ClosedLoopDataLMPC(dt, TimeLMPC, v0,Laps)
LMPCOpenLoopData = LMPCprediction(N, n, d, TimeLMPC, numSS_Points, Laps) #生成所需大小的初始零矩阵
LMPCSimulator = Simulator(map, 1, 1)

LMPController = ControllerLMPC(numSS_Points, numSS_it, N, Qslack, Qlane, Q_LMPC, R_LMPC, dR_LMPC, dt, map, Laps, TimeLMPC, LMPC_Solver, inputConstr)
LMPController.addTrajectory(ClosedLoopDataPID)
LMPController.addTrajectory(ClosedLoopDataLTV_MPC)

x0           = np.zeros((1,n))
x0_glob      = np.zeros((1,n))
x0[0,:]      = ClosedLoopLMPC.x[0,:]
x0_glob[0,:] = ClosedLoopLMPC.x_glob[0,:]

if RunLMPC == 1:
    for it in range(2, Laps):

        ClosedLoopLMPC.updateInitialConditions(x0, x0_glob)
        LMPCSimulator.Sim(ClosedLoopLMPC, LMPController, LMPCOpenLoopData)
        LMPController.addTrajectory(ClosedLoopLMPC)
        ClosedLoopLMPC.SimTimeTot[it]=ClosedLoopLMPC.SimTime
        if LMPController.feasible == 0:
            break
        else:
            # Reset Initial Conditions 重置初始条件
            x0[0,:]      = ClosedLoopLMPC.x[ClosedLoopLMPC.SimTime, :] - np.array([0, 0, 0, 0, map.TrackLength, 0])
            x0_glob[0,:] = ClosedLoopLMPC.x_glob[ClosedLoopLMPC.SimTime, :]

    file_data = open(sys.path[0]+'/data/LMPController.obj', 'wb')
    pickle.dump(ClosedLoopLMPC, file_data)
    pickle.dump(LMPController, file_data)
    pickle.dump(LMPCOpenLoopData, file_data)
    file_data.close()
else:
    file_data = open(sys.path[0]+'/data/LMPController.obj', 'rb')
    ClosedLoopLMPC = pickle.load(file_data)
    LMPController  = pickle.load(file_data)
    LMPCOpenLoopData  = pickle.load(file_data)
    file_data.close()

print "===== LMPC terminated"
# ======================================================================================================================
# ==========================================The Transfer of LMPC  ======================================================
# ======================================================================================================================
#开始迁移学习
print "Starting TLMPC"
ClosedLoopTLMPC = ClosedLoopData(dt, TimeTLMPC, v0)
TLMPCOpenLoopData = TLMPCprediction(N, n, d, TimeTLMPC, numTSS_Points, TLaps,numSR_Points) #生成所需大小的初始零矩阵
#正常仿真
#TLMPCSimulator = Simulator(map, 1, 0,1)
#测试
TLMPCSimulator = testSimulator(map, 1, 0,1)

TLMPController = ControllerTLMPC(numSS_Points, numSR_Points, numSS_it, N, Qslack, Qlane, Q_LMPC, R_LMPC, dR_LMPC, dt, map, Laps, TimeLMPC, LMPC_Solver, inputConstr)
#TLMPController.addTrajectory(ClosedLoopDataPID)
#TLMPController.addTrajectory(ClosedLoopDataLTV_MPC)
TLMPController.addTrajectory(ClosedLoopLMPC,LMPController,Laps)
x0           = np.zeros((1,n))
x0_glob      = np.zeros((1,n))
x0[0,:]      = ClosedLoopTLMPC.x[0,:]
x0_glob[0,:] = ClosedLoopTLMPC.x_glob[0,:]

if RunTLMPC == 1:
    ClosedLoopTLMPC.updateInitialConditions(x0, x0_glob)
    TLMPCSimulator.Sim(ClosedLoopTLMPC, TLMPController, 0,TLMPCOpenLoopData) #闭环数据，控制器，开环预测值和ss集
   # LMPController.addTrajectory(ClosedLoopLMPC)
    file_data = open(sys.path[0]+'/data/TLMPController.obj', 'wb')
    pickle.dump(ClosedLoopTLMPC, file_data)
    pickle.dump(TLMPController, file_data)
    pickle.dump(TLMPCOpenLoopData, file_data)
    file_data.close()
#else:
 #   file_data = open(sys.path[0]+'/data/TLMPController.obj', 'rb')
 #   ClosedLoopTLMPC = pickle.load(file_data)
 #   TLMPController  = pickle.load(file_data)
#    TLMPCOpenLoopData  = pickle.load(file_data)
  #  file_data.close()

print "===== TLMPC terminated"

# ======================================================================================================================
# ========================================= PLOT TRACK =================================================================
# ======================================================================================================================
#for i in range(0, LMPController.it):
#    print "Lap time at iteration ", i, " is ", LMPController.Qfun[0, i]*dt, "s"



print "===== Start Plotting"
if plotFlag == 1:
    plotTrajectory(map, ClosedLoopDataPID.x, ClosedLoopDataPID.x_glob, ClosedLoopDataPID.u)

if plotFlagMPC == 1:
    plotTrajectory(map, ClosedLoopDataLTI_MPC.x, ClosedLoopDataLTI_MPC.x_glob, ClosedLoopDataLTI_MPC.u)

if plotFlagMPC_tv == 1:
    plotTrajectory(map, ClosedLoopDataLTV_MPC.x, ClosedLoopDataLTV_MPC.x_glob, ClosedLoopDataLTV_MPC.u)

if plotFlagLMPC == 1:
    plotClosedLoopLMPC(LMPController, map)

if animation_xyFlag == 1:
    #animation_xy(map, LMPCOpenLoopData, LMPController, 3)
    saveGif_xyResults(map, LMPCOpenLoopData, LMPController,40 )
    plt.show()
if animation_stateFlag == 1:
   # animation_states(map, LMPCOpenLoopData, LMPController, 5)
    Save_statesAnimation(map, LMPCOpenLoopData, LMPController, 5)
    plt.show()
# unityTestChangeOfCoordinates(map, ClosedLoopDataPID)
# unityTestChangeOfCoordinates(map, ClosedLoopDataLTI_MPC)
# unityTestChangeOfCoordinates(map, ClosedLoopLMPC)
if plotFlagTLMPC==1:
    TLMPC_xyResults(map,TLMPCOpenLoopData,TLMPController,0)
    plt.show()


