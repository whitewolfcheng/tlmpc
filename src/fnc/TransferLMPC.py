#coding=utf-8

import pdb
import numpy as np
from cvxopt import spmatrix, matrix, solvers
from numpy import linalg as la
from scipy import linalg
from scipy import sparse
from cvxopt.solvers import qp
import datetime
from Utilities import Curvature
from numpy import hstack, inf, ones
from scipy.sparse import vstack

class ControllerTLMPC():
    """Create the LMPC
    Attributes:
        solve: given x0 computes the control action
        addTrajectory: given a ClosedLoopData object adds the trajectory to SS, Qfun, uSS and updates the iteration index
        addPoint: this function allows to add the closed loop data at iteration j to the SS of iteration (j-1)
        update: this function can be used to set SS, Qfun, uSS and the iteration index.
    """

    def __init__(self, numSS_Points,numSR_Points, numSS_it, N, Qslack, Qlane, Q, R, dR, dt,  map, Laps, TimeLMPC, Solver, inputConstraints):
        """Initialization
        Arguments:
            numSS_Points: number of points selected from the previous trajectories to build SS
            numSS_it: number of previois trajectories selected to build SS
            N: horizon length
            Q,R: weight to define cost function h(x,u) = ||x||_Q + ||u||_R
            dR: weight to define the input rate cost h(x,u) = ||x_{k+1}-x_k||_dR
            n,d: state and input dimensiton
            shift: given the closest point x_t^j to x(t) the controller start selecting the point for SS from x_{t+shift}^j
                  #x_{t+shift}^j   第j次迭代t+shift时刻离现在x（t）的点最近
            map: map
            Laps: maximum number of laps the controller can run (used to avoid dynamic allocation)（避免动态分配……应该是动态内存分配之类的）
            TimeLMPC: maximum time [s] that an lap can last (used to avoid dynamic allocation)
            Solver: solver used in the reformulation of the LMPC as QP
        """
        self.numSS_Points = numSS_Points
        self.numSR_points=numSR_Points
        self.numSS_it     = numSS_it
        self.N = N
        self.Qslack = Qslack
        self.Qlane = Qlane
        self.Q = Q
        self.R = R
        self.dR = dR
        self.n = Q.shape[1]
        self.d = R.shape[1]
        self.dt = dt
        self.map = map
        self.Solver = Solver
        self.LapTime = 0
        self.itUsedSysID = 1
        self.inputConstraints = inputConstraints

        self.OldInput = np.zeros((1,2))
#初始化
        # Initialize the following quantities to avoid dynamic allocation
        NumPoints = int(TimeLMPC / dt) + 1
        self.LapCounter = 10000 * np.ones(Laps).astype(int)        # Time at which each j-th iteration is completed
        self.TimeSS     = 10000 * np.ones(Laps).astype(int)        # Time at which each j-th iteration is completed
        self.SS         = 10000 * np.ones((NumPoints, 6, Laps))    # Sampled Safe SS
        self.uSS        = 10000 * np.ones((NumPoints, 2, Laps))    # Input associated with the points in SS
        self.bcur       = 100000 * np.ones((NumPoints, Laps))       # bcur: 数据集中各点处的曲率
        self.arclen = 0    * np.ones((NumPoints, Laps))            # arclen:该点与前一点之间的弧长
        self.arclens = 0 * np.ones((NumPoints, Laps))               # arclens:起点到该点间的弧长
        self.SS_glob    = 10000 * np.ones((NumPoints, 6, Laps))    # SS in global (X-Y) used for plotting

        self.zVector = np.array([0.0, 0.0, 0.0, 0.0, 10.0, 0.0])

        #初始化迭代次数
        # Initialize the controller iteration
        self.it      = 0
        #LMPC不等式约束矩阵
        # Build matrices for inequality constraints
    #    self.F, self.b = _LMPC_BuildMatIneqConst(self)

        self.xPred = []

    def solve(self, x0, uOld=np.zeros([0, 0])):
        """Computes control action
        Arguments:
            x0: current state position
        """
        n = self.n;
        d = self.d
#        F = self.F;
 #       b = self.b
        SS = self.SS;
        bcur = self.bcur
        arclen=self.arclen
        arclens=self.arclens
        uSS = self.uSS;
        TimeSS = self.TimeSS
        Q = self.Q;
        R = self.R
        dR = self.dR;
        OldInput = self.OldInput
        N = self.N;
        dt = self.dt
        it = self.it
        numSS_Points = self.numSS_Points
        numSR_Points=self.numSR_points
        Qslack = self.Qslack

        LinPoints = self.LinPoints
        LinInput = self.LinInput
        map = self.map
        PointAndTangent = map.PointAndTangent

        # Select Points from SS
        if (self.zVector[4] - x0[4] > map.TrackLength / 2):
            self.zVector[4] = np.max([self.zVector[4] - map.TrackLength, 0])
            self.LinPoints[4, -1] = self.LinPoints[4, -1] - map.TrackLength
        #判断相似路径
        SR_PointSelectedTot = np.empty((n, 0))
        SS_PointSelectedTot = np.empty((n, 0))
        Succ_SS_PointSelectedTot = np.empty((n, 0))
        Succ_uSS_PointSelectedTot = np.empty((d, 0))
        Qfun_SelectedTot = np.empty((0))

        maxsspoint = np.argmax(TimeSS[0])
        Mindiff=100.0
        for ii in range(0,maxsspoint-numSR_Points):
            #记忆中的numSR_Points个点的曲率
            curb=bcur[ii:ii+numSR_Points,:]
            #当前路径中心线的numSR_Points个点的对应s值及曲率
            nowarclen=arclens[ii:ii+numSR_Points,:]
            nowxs=self.zVector[4]  #修正s值
            diffarclens=nowxs-nowarclen[9,:]
            nowcarclens=nowarclen+diffarclens
            nowcur = np.empty((nowarclen.shape[0], nowarclen[1]))
            for i in range(0,nowarclen.shape[0]):
                for j in range(0, nowarclen.shape[1]):
                    nowcur[i,j] = Curvature(nowcarclens[i,j], PointAndTangent)
            diffcur=nowcur-curb
            diffcurnorm = la.norm(diffcur, 1, axis=0)
            MinJ = np.argmin(diffcurnorm)
            if Mindiff>=diffcurnorm[MinJ]:
                MinJJ=MinJ
                MinII=ii
                Mindiff=diffcurnorm[MinJ]
        similarroadpoint=SS[MinII:MinII+numSR_Points,:, MinJJ]  #最相似路段的轨迹点
        self.TSR_PointSelectedTot = similarroadpoint
            







    def addTrajectory(self, ClosedLoopData):
        """update iteration index and construct SS, uSS and Qfun
        Arguments:
            ClosedLoopData: ClosedLoopData object
        """
        it = self.it

        self.TimeSS[it] = ClosedLoopData.SimTime
        self.LapCounter[it] = ClosedLoopData.SimTime
        self.SS[0:(self.TimeSS[it] + 1), :, it] = ClosedLoopData.x[0:(self.TimeSS[it] + 1), :]
        self.SS_glob[0:(self.TimeSS[it] + 1), :, it] = ClosedLoopData.x_glob[0:(self.TimeSS[it] + 1), :]
        self.uSS[0:self.TimeSS[it], :, it] = ClosedLoopData.u[0:(self.TimeSS[it]), :]

        for i in np.arange(1, self.SS_glob.shape[0]-2):
            x = self.SS_glob[i-1:i + 2, 4, it]
            y = self.SS_glob[i-1:i + 2, 5, it]
            cur=PJcurvature(x,y)
            #求曲率
            self.bcur[i, it] = cur
            #求弧长
            bcur1=cur
            if i==1:
                self.bcur[0,it]=self.bcur[1,it]
                bcur0=bcur1
            else :
                bcur0=self.bcur[i-1,it]
            arclen=Arclength(bcur1,bcur0,x[1],y[1],x[0],y[0])
            self.arclen[i,it]=arclen
            self.arclens[i,it]=self.arclens[i-1,it]+self.arclen[i,it]

       # self.bcur[-1,it]=self.bcur[-2,it]
        if self.it == 0:
            self.LinPoints = self.SS[1:self.N + 2, :, it]
            self.LinInput = self.uSS[1:self.N + 1, :, it]

        self.it = self.it + 1

    def addPoint(self, x, u,xglob):
        """at iteration j add the current point to SS, uSS and Qfun of the previous iteration
        Arguments:
            x: current state
            u: current input
            i: at the j-th iteration i is the time at which (x,u) are recorded
        """
        Counter = self.TimeSS[self.it - 1]
        self.SS[Counter, :, self.it - 1] = x + np.array([0, 0, 0, 0, self.map.TrackLength, 0])
        self.uSS[Counter, :, self.it - 1] = u
        self.SS_glob[Counter, :, self.it - 1]=xglob
        x = self.SS_glob[Counter-2:Counter+1, 4, self.it-1]
        y = self.SS_glob[Counter-2:Counter+1, 5, self.it-1]
        cur=PJcurvature(x,y)
        self.bcur[Counter-1,self.it-1]=cur
        arclen = Arclength(cur, self.bcur[Counter-2,self.it-1], x[1], y[1], x[0], y[0])
        self.arclen[Counter-1, self.it-1] = arclen
        self.arclens[Counter-1, self.it-1] = self.arclens[Counter-2, self.it-1] + self.arclen[Counter-1, self.it-1]
        # The above  are needed as the once the predicted trajectory has crossed the finish line the goal is
        # to reach the end of the lap which is about to start

        self.TimeSS[self.it - 1] = self.TimeSS[self.it - 1] + 1

#求相邻三点的中间点的曲率
def PJcurvature(x, y):
    """
    input  : the coordinate of the three point
    output : the curvature and norm direction
    refer to https://github.com/Pjer-zhang/PJCurvature for detail
    """
    #如果两点共线
    if abs((y[2]-y[0])*(x[1]-x[0])-(y[1]-y[0])*(x[2]-x[0]))<=1e-6:
        kappa=0
    else:
        t_a = la.norm([x[1] - x[0], y[1] - y[0]])
        t_b = la.norm([x[2] - x[1], y[2] - y[1]])

        M = np.array([
            [1, -t_a, t_a ** 2],
            [1, 0, 0],
            [1, t_b, t_b ** 2]
            ])

        a = np.matmul(la.inv(M), x)
        b = np.matmul(la.inv(M), y)

        kappa = 2 * (a[2] * b[1] - b[2] * a[1]) / (a[1] ** 2. + b[1] ** 2.) ** (1.5)
    return kappa
#求两点间弧长
def Arclength(bcur1,bcur0,x1,y1,x0,y0):
    bcur=abs(0.5*(bcur1+bcur0))
    if bcur<=1e-6:
        arclen=np.sqrt((x1-x0)**2+(y1-y0)**2)
    else:
        r=1/bcur
        d=np.sqrt((x1-x0)**2+(y1-y0)**2)
        rad=np.arcsin(0.5*d*bcur)
        arclen=rad*r
    return arclen
