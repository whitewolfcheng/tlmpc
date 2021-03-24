# coding=utf-8

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

# from osqp import OSQP

solvers.options['show_progress'] = False


class NewControllerLMPC():
    """Create the LMPC
    Attributes:
        solve: given x0 computes the control action
        addTrajectory: given a ClosedLoopData object adds the trajectory to SS, Qfun, uSS and updates the iteration index
        addPoint: this function allows to add the closed loop data at iteration j to the SS of iteration (j-1)
        update: this function can be used to set SS, Qfun, uSS and the iteration index.
    """

    def __init__(self, numSS_Points, numSS_it, N, Qslack, Qlane, Q, R, dR, dt, map, Laps, TLaps,TimeLMPC, Solver,
                 inputConstraints):
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
        self.numSS_it = numSS_it
        self.TLaps=TLaps
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
        self.itUsedSysID = 2
        self.inputConstraints = inputConstraints

        self.OldInput = np.zeros((1, 2))
        # 初始化
        # Initialize the following quantities to avoid dynamic allocation
        NumPoints = int(TimeLMPC / dt) + 1
        self.transSSTime  = 0*np.ones(NumPoints).astype(int)  #迁移后的SS集中的点个数(分段统计)
        self.tSSlapsTime= 0*np.ones(TLaps).astype(int)  #迁移后的SS集中每一圈的点个数
        self.LapCounter = 10000 * np.ones(Laps).astype(int)  # Time at which each j-th iteration is completed 只到终点
        self.TimeSS = 10000 * np.ones(Laps).astype(int)  # Time at which each j-th iteration is completed  到终点后继续增加下一圈的点
        self.TimeTSS=10000 * np.ones(Laps).astype(int) #存放在跑新的一圈之后用正常LMPC算法存放的数据
        self.SS = 10000 * np.ones((NumPoints, 6, Laps))  # Sampled Safe SS 用来建立跑第一圈所需要的SS集，即存放迁移后的数据
        self.tSS = 10000 * np.ones((NumPoints, 6, Laps))  # 用来储存第一圈实际跑的数据
        self.uSS = 10000 * np.ones((NumPoints, 2, Laps))  # Input associated with the points in SS
        self.tSSparaTime = 0 * np.ones((NumPoints, TLaps))  # 迁移后，每一圈里各片段的点个数
        self.tuSS=10000 * np.ones((NumPoints, 2, Laps))  # Input associated with the points in SS
        self.Qfun = 0 * np.ones((NumPoints, Laps))  # Qfun: cost-to-go from each point in SS
        self.SS_glob = 10000 * np.ones((NumPoints, 6, Laps))  # SS in global (X-Y) used for plotting
        self.zVector = np.array([0.0, 0.0, 0.0, 0.0, 10.0, 0.0])

        # 初始化迭代次数
        # Initialize the controller iteration
        self.it = 0
        # LMPC不等式约束矩阵
        # Build matrices for inequality constraints
        self.F, self.b = _LMPC_BuildMatIneqConst(self)

        self.xPred = []

    def solve(self, x0, uOld=np.zeros([0, 0])):
        """Computes control action
        Arguments:
            x0: current state position
        """
        n = self.n;
        d = self.d
        F = self.F;
        b = self.b
        SS = self.SS;
        Qfun = self.Qfun
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
        Qslack = self.Qslack

        LinPoints = self.LinPoints
        LinInput = self.LinInput   
        map = self.map
        # Select Points from SS
        if (self.zVector[4] - x0[4] > map.TrackLength / 2):
            self.zVector[4] = np.max([self.zVector[4] - map.TrackLength, 0])
            self.LinPoints[4, -1] = self.LinPoints[4, -1] - map.TrackLength


        # argsort函数返回的是数组值从小到大的索引值， Qfun为ss中从某点开始的代价函数
    #    sortedLapTime = np.argsort(self.Qfun[0, 0:it])

        SS_PointSelectedTot = np.empty((n, 0))
        Succ_SS_PointSelectedTot = np.empty((n, 0))
        Succ_uSS_PointSelectedTot = np.empty((d, 0))
        Qfun_SelectedTot = np.empty((0))
        SS_PointSelected=np.empty((n,numSS_Points / self.numSS_it+1,self.TLaps))
        uSS_PointSelected=np.empty((d,numSS_Points / self.numSS_it+1,self.TLaps))
        Qfun_Selected=np.empty((numSS_Points / self.numSS_it+1,self.TLaps))
        MinX0diff=np.empty((self.TLaps))
        # 将符合条件的点加 到ss集和qfun代价函数中
        # print(sortedLapTime[0:self.numSS_it])
        for jj in range(0,self.TLaps): #先遍历所有圈数挑选每一圈的SS集
            # print(jj)
            SS_PointSelected[:,:,jj], uSS_PointSelected[:,:,jj], Qfun_Selected[:,jj], MinX0diff[jj] = _SelectPoints(self, jj, self.zVector,
                                                                               numSS_Points / self.numSS_it + 1) #Mindiff 是选中点与当前点的距离
            # print(SS_PointSelected.shape)
        #    print 'MinX0diff' ,MinX0diff[jj]
        sortedLapTime=np.argsort(MinX0diff)  #根据和当前点的距离对挑选的SS集进行从近到远的排序
     #   with open('Mindiff.csv', 'ab') as f:
  #          np.savetxt(f, MinX0diff[sortedLapTime[0:2]], fmt='%f')
        for jjs in sortedLapTime[0:self.numSS_it]: #根据所需SS集段数的要求，挑选最近的SS集
            # 下一状态
        #    print '选中迭代圈数',jjs
            Succ_SS_PointSelectedTot = np.append(Succ_SS_PointSelectedTot, SS_PointSelected[:, 1:,jjs], axis=1)
            Succ_uSS_PointSelectedTot = np.append(Succ_uSS_PointSelectedTot, uSS_PointSelected[:, 1:,jjs], axis=1)
            # 当前状态
            # print(Succ_SS_PointSelectedTot)
            SS_PointSelectedTot = np.append(SS_PointSelectedTot, SS_PointSelected[:, 0:-1,jjs], axis=1)
            Qfun_SelectedTot = np.append(Qfun_SelectedTot, Qfun_Selected[0:-1,jjs], axis=0)

        self.SS_PointSelectedTot = SS_PointSelectedTot
        self.Qfun_SelectedTot = Qfun_SelectedTot
        # print(self.zVector[4])
        # print(SS_PointSelectedTot[4,0:20])
        #  print(SS_PointSelectedTot[4, 20:40])
        # 等式约束的建立
        startTimer = datetime.datetime.now()
        Atv, Btv, Ctv, indexUsed_list = _LMPC_EstimateABC(self, sortedLapTime)
        endTimer = datetime.datetime.now();
        deltaTimer = endTimer - startTimer
        L, npG, npE = _LMPC_BuildMatEqConst(self, Atv, Btv, Ctv, N, n, d)
        self.linearizationTime = deltaTimer
        # 建立终端代价约束
        # Build Terminal cost and Constraint
        G, E = _LMPC_TermConstr(self, npG, npE, N, n, d, SS_PointSelectedTot)
        M, q = _LMPC_BuildMatCost(self, Qfun_SelectedTot, numSS_Points, N, Qslack, Q, R, dR, OldInput)

        # Solve QP 两种求解方法：CVX法和OSQP法
        startTimer = datetime.datetime.now()

        if self.Solver == "CVX":
            res_cons = qp(M, matrix(q), F, matrix(b), G, E * matrix(x0) + L)
            #调查等式约束
         #   np.savetxt('G.csv', G, delimiter=',', fmt='%f')
         #   np.savetxt('E.csv', E * matrix(x0) + L, delimiter=',', fmt='%f')
            if res_cons['status'] == 'optimal':
                feasible = 1
            else:
                feasible = 0
                print res_cons['status']
            Solution = np.squeeze(res_cons['x'])

        elif self.Solver == "OSQP":
            # Adaptation for QSQP from https://github.com/alexbuyval/RacingLMPC/
            res_cons, feasible = osqp_solve_qp(sparse.csr_matrix(M), q, sparse.csr_matrix(F), b, sparse.csr_matrix(G),
                                               np.add(np.dot(E, x0), L[:, 0]))
            Solution = res_cons.x

        self.feasible = feasible

        endTimer = datetime.datetime.now();
        deltaTimer = endTimer - startTimer
        self.solverTime = deltaTimer

        # Extract solution and set linerizations points
        xPred, uPred, lambd, slack = _LMPC_GetPred(Solution, n, d, N, np)
  #      with open('zhongduan.csv', 'ab') as f:
 #           np.savetxt(f, np.dot(SS_PointSelectedTot, lambd), fmt='%f')
  #      with open('zVector.csv', 'ab') as f:
   #         np.savetxt(f, self.zVector, fmt='%f')
    #    with open('diffzhongduan.csv', 'ab') as f:
    #        np.savetxt(f, np.dot(SS_PointSelectedTot, lambd)-self.zVector, fmt='%f')
        self.zVector = np.dot(Succ_SS_PointSelectedTot, lambd)
        self.uVector = np.dot(Succ_uSS_PointSelectedTot, lambd)
        self.xPred = xPred.T
      #  print self.zVector,self.uVector
        if self.N == 1:
            self.uPred = np.array([[uPred[0], uPred[1]]])
            self.LinInput = np.array([[uPred[0], uPred[1]]])
        else:
            self.uPred = uPred.T
            self.LinInput = np.vstack((uPred.T[1:, :], self.uVector))

        self.LinPoints = np.vstack((xPred.T[1:, :], self.zVector))
        print '状态',self.LinPoints[:,4]
   #     print '输入',self.LinInput

        self.OldInput = uPred.T[0, :]

    def firstTrajectory(self, ClosedLoopData):
        """update iteration index and construct SS, uSS and      Qfun
        Arguments:
            ClosedLoopData: ClosedLoopData object
        """
        TLaps = ClosedLoopData.TLaps
        self.tSSparaTime=ClosedLoopData.tSSparaTime
        for i in range(0,TLaps):
            if i == 0:
                tSSlapslast=0
            else:
                tSSlapslast=ClosedLoopData.tSSlapsTime[i-1]
            self.TimeSS[i]=ClosedLoopData.tSSlapsTime[i]-tSSlapslast  #每一圈有多少个点
            self.SS[0:self.TimeSS[i], :, i] = ClosedLoopData.transSS[tSSlapslast:ClosedLoopData.tSSlapsTime[i], :]
            self.uSS[0:self.TimeSS[i], :, i] = ClosedLoopData.transuSS[tSSlapslast:ClosedLoopData.tSSlapsTime[i], :]



        if self.it == 0:
            self.LinPoints = self.SS[1:self.N + 2, :, 0]
            self.LinInput = self.uSS[1:self.N + 1, :, 0]
           # print self.LinPoints[:,4]
      #      self.zVector =self.SS[self.N,:,0]




    def addPoint(self, x, u, xglob):
        """at iteration j add the current point to SS, uSS and Qfun of the previous iteration
        Arguments:
            x: current state
            u: current inputsd
            i: at the j-th iteration i is the time at which (x,u) are recorded
        """
        Counter = self.TimeSS[self.it - 1]
        print Counter
        print self.it
        print self.TimeSS
        self.SS[Counter, :, self.it - 1] = x + np.array([0, 0, 0, 0, self.map.TrackLength, 0])
        self.uSS[Counter, :, self.it - 1] = u
        self.SS_glob[Counter, :, self.it - 1] = xglob

        # The above two lines are needed as the once the predicted trajectory has crossed the finish line the goal is
        # to reach the end of the lap which is about to start
        if self.Qfun[Counter, self.it - 1] == 0:
            self.Qfun[Counter, self.it - 1] = self.Qfun[Counter, self.it - 1] - 1

        self.TimeSS[self.it - 1] = self.TimeSS[self.it - 1] + 1

    def update(self, SS, uSS, Qfun, TimeSS, it, LinPoints, LinInput):
        """update controller parameters. This function is useful to transfer information among LMPC controller
           with different tuning
        Arguments:
            SS: sampled safe set
            uSS: input associated with the points in SS
            Qfun: Qfun: cost-to-go from each point in SS
            TimeSS: time at which each j-th iteration is completed
            it: current iteration
            LinPoints: points used in the linearization and system identification procedure
            LinInput: inputs associated with the points used in the linearization and system identification procedure
        """
        self.SS = SS
        self.uSS = uSS
        self.Qfun = Qfun
        self.TimeSS = TimeSS
        self.it = it

        self.LinPoints = LinPoints
        self.LinInput = LinInput


# ======================================================================================================================
# ======================================================================================================================
# =============================== Internal functions for LMPC reformulation to QP ======================================
# ======================================================================================================================
# ======================================================================================================================

# 代价函数
def _LMPC_BuildMatCost(LMPC, Sel_Qfun, numSS_Points, N, Qslack, Q, R, dR, uOld):
    n = Q.shape[0]
    P = Q
    vt = 2
    Qlane = LMPC.Qlane

    b = [Q] * (N)
    Mx = linalg.block_diag(*b)

    c = [R + 2 * np.diag(dR)] * (N)  # Need to add dR for the derivative input cost

    Mu = linalg.block_diag(*c)
    # Need to condider that the last input appears just once in the difference
    Mu[Mu.shape[0] - 1, Mu.shape[1] - 1] = Mu[Mu.shape[0] - 1, Mu.shape[1] - 1] - dR[1]
    Mu[Mu.shape[0] - 2, Mu.shape[1] - 2] = Mu[Mu.shape[0] - 2, Mu.shape[1] - 2] - dR[0]

    # Derivative Input Cost
    OffDiaf = -np.tile(dR, N - 1)  # 将dR横着复制N-1次
    np.fill_diagonal(Mu[2:], OffDiaf)
    np.fill_diagonal(Mu[:, 2:], OffDiaf)
    # np.savetxt('Mu.csv', Mu, delimiter=',', fmt='%f')

    M00 = linalg.block_diag(Mx, P, Mu)
    quadLaneSlack = Qlane[0] * np.eye(2 * LMPC.N)
    M0 = linalg.block_diag(M00, np.zeros((numSS_Points, numSS_Points)), Qslack, quadLaneSlack)
    # np.savetxt('M0.csv', M0, delimiter=',', fmt='%f')

    xtrack = np.array([vt, 0, 0, 0, 0, 0])
    q0 = - 2 * np.dot(np.append(np.tile(xtrack, N + 1), np.zeros(R.shape[0] * N)), M00)

    # Derivative Input
    q0[n * (N + 1):n * (N + 1) + 2] = -2 * np.dot(uOld, np.diag(dR))

    # np.savetxt('q0.csv', q0, delimiter=',', fmt='%f')
    linLaneSlack = Qlane[1] * np.ones(2 * LMPC.N)

    q = np.append(np.append(np.append(q0, Sel_Qfun), np.zeros(Q.shape[0])), linLaneSlack)

    # np.savetxt('q.csv', q, delimiter=',', fmt='%f')

    M = 2 * M0  # Need to multiply by two because CVX considers 1/2 in front of quadratic cost

    if LMPC.Solver == "CVX":
        M_sparse = spmatrix(M[np.nonzero(M)], np.nonzero(M)[0].astype(int), np.nonzero(M)[1].astype(int), M.shape)
        M_return = M_sparse
    else:
        M_return = M

    return M_return, q


# 不等式约束
def _LMPC_BuildMatIneqConst(LMPC):
    N = LMPC.N
    n = LMPC.n
    numSS_Points = LMPC.numSS_Points
    # Buil the matrices for the state constraint in each region. In the region i we want Fx[i]x <= bx[b]
    Fx = np.array([[0., 0., 0., 0., 0., 1.],
                   [0., 0., 0., 0., 0., -1.]])

    bx = np.array([[LMPC.map.halfWidth],  # max ey
                   [LMPC.map.halfWidth]])  # max ey

    # Buil the matrices for the input constraint in each region. In the region i we want Fx[i]x <= bx[b]
    Fu = np.array([[1., 0.],
                   [-1., 0.],
                   [0., 1.],
                   [0., -1.]])

    bu = np.array([[LMPC.inputConstraints[0, 0]],  # Max Steering
                   [LMPC.inputConstraints[0, 1]],  # Min Steering
                   [LMPC.inputConstraints[1, 0]],  # Max Acceleration
                   [LMPC.inputConstraints[1, 1]]])  # Min Acceleration

    # Now stuck the constraint matrices to express them in the form Fz<=b. Note that z collects states and inputs
    # Let's start by computing the submatrix of F relates with the state
    rep_a = [Fx] * (N)
    Mat = linalg.block_diag(*rep_a)
    NoTerminalConstr = np.zeros((np.shape(Mat)[0], n))  # No need to constraint also the terminal point
    Fxtot = np.hstack((Mat, NoTerminalConstr))
    bxtot = np.tile(np.squeeze(bx), N)

    # Let's start by computing the submatrix of F relates with the input
    rep_b = [Fu] * (N)
    Futot = linalg.block_diag(*rep_b)
    butot = np.tile(np.squeeze(bu), N)

    # Let's stack all together
    rFxtot, cFxtot = np.shape(Fxtot)
    rFutot, cFutot = np.shape(Futot)
    # 横向拼接
    Dummy1 = np.hstack((Fxtot, np.zeros((rFxtot, cFutot))))
    Dummy2 = np.hstack((np.zeros((rFutot, cFxtot)), Futot))

    FDummy = np.vstack((Dummy1, Dummy2))
    I = -np.eye(numSS_Points)
    FDummy2 = linalg.block_diag(FDummy, I)
    Fslack = np.zeros((FDummy2.shape[0], n))
    F_hard = np.hstack((FDummy2, Fslack))

    LaneSlack = np.zeros((F_hard.shape[0], 2 * N))
    colIndexPositive = []
    rowIndexPositive = []
    colIndexNegative = []
    rowIndexNegative = []
    for i in range(0, N):
        colIndexPositive.append(i * 2 + 0)
        colIndexNegative.append(i * 2 + 1)

        rowIndexPositive.append(i * Fx.shape[0] + 0)  # Slack on second element of Fx
        rowIndexNegative.append(i * Fx.shape[0] + 1)  # Slack on third element of Fx

    LaneSlack[rowIndexPositive, colIndexPositive] = -1.0
    LaneSlack[rowIndexNegative, rowIndexNegative] = -1.0

    F_1 = np.hstack((F_hard, LaneSlack))

    I = - np.eye(2 * N)
    Zeros = np.zeros((2 * N, F_hard.shape[1]))
    Positivity = np.hstack((Zeros, I))

    F = np.vstack((F_1, Positivity))

    # np.savetxt('F.csv', F, delimiter=',', fmt='%f')
    # pdb.set_trace()

    b_1 = np.hstack((bxtot, butot, np.zeros(numSS_Points)))

    b = np.hstack((b_1, np.zeros(2 * N)))
    # np.savetxt('b.csv', b, delimiter=',', fmt='%f')

    if LMPC.Solver == "CVX":
        F_sparse = spmatrix(F[np.nonzero(F)], np.nonzero(F)[0].astype(int), np.nonzero(F)[1].astype(int), F.shape)
        F_return = F_sparse
    else:
        F_return = F

    return F_return, b


def _SelectPoints(LMPC, it, x0, numSS_Points):
    SS = LMPC.SS
    uSS = LMPC.uSS
    TimeSS = LMPC.TimeSS
    tSSparaTime=LMPC.tSSparaTime
    Qfun = LMPC.Qfun
    xPred = LMPC.xPred
    map = LMPC.map
    TrackLength = map.TrackLength
    currIt = LMPC.it
    x = SS[0:TimeSS[it] , :, it]
    u = uSS[0:TimeSS[it] , :, it]
    oneVec = np.ones((x.shape[0], 1))
    x0Vec = (np.dot(np.array([x0]).T, oneVec.T)).T  # 每行是一个x0，共x.shape[0]行
    diff = x - x0Vec
    norm = la.norm(diff, 1, axis=1)
    MinNorm = np.argmin(norm)
#    Mindiff=norm[MinNorm]
    Mindiff=la.norm(x[MinNorm,:]-x0,1)
    #找到离当前终端约束点最近的点在哪个迁移片段内
    paradiff=np.squeeze(MinNorm-tSSparaTime[:,it])
    indexpara1=np.where(paradiff<0)[0][0]  #右边的点
    if indexpara1==0:
        para1=0
    else :
        para1=int(tSSparaTime[indexpara1-1,it]) #左边界
    para2=int(tSSparaTime[indexpara1,it]-1)   #右边界
    # 如果当前点也在下一段轨迹中，并且还超了三个点,提前跳到下一路段
    if para2+1+5<x.shape[0]:
      if x[para2+1+5,4]<=x0[4] :
        para1=para2+1
        para2=int(tSSparaTime[indexpara1+1,it]-1)
        xlocal=SS[para1:para2+1 , :, it]
        oneVeclocal=np.ones((para2-para1+1,1))
        x0Vec=(np.dot(np.array([x0]).T, oneVeclocal.T)).T
        difflocal=xlocal-x0Vec
        normlocal=la.norm(difflocal,1,axis=1)
        MinNorm=np.argmin(normlocal)+para1
        Mindiff = la.norm(x[MinNorm, :] - x0, 1)
 #   if MinNorm==para2: #准备跳车
     #   print "开始跳车"
  #      MinNorm = MinNorm+1  #往前跨一步
    #如果最临近点在该圈的最后一个点
    if MinNorm>=x.shape[0] and it<LMPC.TLaps-1:
            MinNorm=0
            x = SS[0:TimeSS[it+1], :, it+1]
            u = uSS[0:TimeSS[it+1], :, it+1]
            x[:,4]=x[:,4]+map.TrackLength
            para1=0
            para2 = int(tSSparaTime[0, it + 1]-1)
            Mindiff=la.norm(x[0,:]-x0,1)
    if MinNorm >= x.shape[0] and it == LMPC.TLaps - 1:
            MinNorm = 0
            x = SS[0:TimeSS[it], :, it ]
            u = uSS[0:TimeSS[it], :, it ]
            x[:, 4] = x[:, 4] + map.TrackLength
            para1 = 0
            para2 = int(tSSparaTime[0, it] - 1)
            Mindiff = la.norm(x[0, :] - x0,1)
 #       else:
  #          para1=int(tSSparaTime[indexpara1,it] )#下一段路径的左边界
   #         para2 =int( tSSparaTime[indexpara1 + 1, it] - 1 )#下一段路径的右边界
  #          Mindiff=norm[MinNorm]
   # print Mindiff


    #从迁移路段中挑选SS集
    #降低SS集断层较大的路段被选中的优先级
    if (MinNorm-para1<numSS_Points/2): #and para1!=0:
        indexSSandQfun=range(para1, para1 + numSS_Points) #如果当前点很靠近迁移路段的最左边，则直接从最左边开始选SS集
   #    Mindiff=Mindiff+la.norm(x[para1,:]-x[para1-1,:],1)
       # print '左边跳车',it
    if (para2-MinNorm<numSS_Points/2):
     #   print '右边跳车',it
        indexSSandQfun=range(para2-numSS_Points+1,para2+1) #如果当前点很靠近迁移路段最右边，则直接从最右边开始选SS集。
    #   Mindiff = Mindiff + la.norm(x[para2, :] - x[para2 + 1, :], 1)
        # 本来python左开右闭原则要两端+1，由于要留出下一步的SS集因此还得两端-1，因此抵消，不用+1也不用-1
    if (MinNorm-para1>=numSS_Points/2)&(para2-MinNorm>=numSS_Points/2):
     #   print '爷不跳了',it
        indexSSandQfun=range(-numSS_Points / 2 + MinNorm, numSS_Points / 2 + MinNorm) #如果当前点在中间，则两边各取一半建立SS集

    #选择SS集
   # if (MinNorm - numSS_Points/2 >= 0):
   #     indexSSandQfun = range(-numSS_Points/2 + MinNorm, numSS_Points/2 + MinNorm)
  # else:
  #      indexSSandQfun = range(MinNorm, MinNorm + numSS_Points)
    SS_Points = x[indexSSandQfun, :].T
    # print(SS_Points.shape)
    SSu_Points = u[indexSSandQfun, :].T
    #   print(SSu_Points.shape)
    Sel_Qfun = (x0[4]-x[indexSSandQfun, 4] ) *10#用SS集中各点S值领先上一个终端约束点的程度代替Q值
  #  print Sel_Qfun,x0[4],x[indexSSandQfun, 4],MinNorm,indexSSandQfun
  #  print Sel_Qfun[Sel_Qfun>0]

    return SS_Points, SSu_Points, Sel_Qfun, Mindiff


# 没有到最终点，则代价+1
def _ComputeCost(x, u, TrackLength):
    Cost = 10000 * np.ones((x.shape[0]))  # The cost has the same elements of the vector x --> time +1
    # Now compute the cost moving backwards in a Dynamic Programming (DP) fashion.
    # We start from the last element of the vector x and we sum the running cost
    for i in range(0, x.shape[0]):
        if (i == 0):  # Note that for i = 0 -- > pick the latest element of the vector x
            Cost[x.shape[0] - 1 - i] = 0
        elif x[x.shape[0] - 1 - i, 4] < TrackLength:
            Cost[x.shape[0] - 1 - i] = Cost[x.shape[0] - 1 - i + 1] + 1
        else:
            Cost[x.shape[0] - 1 - i] = 0

    return Cost


def _LMPC_TermConstr(LMPC, G, E, N, n, d, SS_Points):
    # 终端约束，让每次滚动的最后一步都回到ss集
    # Update the matrices for the Equality constraint in the LMPC. Now we need an extra row to constraint the terminal point to be equal to a point in SS
    # The equality constraint has now the form: G_LMPC*z = E_LMPC*x0 + TermPoint.
    # Note that the vector TermPoint is updated to constraint the predicted trajectory into a point in SS. This is done in the FTOCP_LMPC function

    TermCons = np.zeros((n, (N + 1) * n + N * d))
    TermCons[:, N * n:(N + 1) * n] = np.eye(n)  # 第N+1步的状态变量系数为1

    G_enlarged = np.vstack((G, TermCons))

    G_lambda = np.zeros((G_enlarged.shape[0], SS_Points.shape[1] + n))
    G_lambda[G_enlarged.shape[0] - n:G_enlarged.shape[0], :] = np.hstack((-SS_Points, np.eye(n)))

    G_LMPC0 = np.hstack((G_enlarged, G_lambda))  # 横向连接
    G_ConHull = np.zeros((1, G_LMPC0.shape[1]))
    G_ConHull[-1, G_ConHull.shape[1] - SS_Points.shape[1] - n:G_ConHull.shape[1] - n] = np.ones((1, SS_Points.shape[1]))
    # 等式左边ss集中每个点的权重，权重之和等于等式右边的L矩阵的最后一个元素，即等于1，# 一般是某个点的权重为1，其余为0，为1的那个点即是预测时域最后一步回到的ss集中的那个点

    G_LMPC_hard = np.vstack((G_LMPC0, G_ConHull))

    SlackLane = np.zeros((G_LMPC_hard.shape[0], 2 * N))
    G_LMPC = np.hstack((G_LMPC_hard, SlackLane))

    E_LMPC = np.vstack((E, np.zeros((n + 1, n))))

    # np.savetxt('G.csv', G_LMPC, delimiter=',', fmt='%f')
    # np.savetxt('E.csv', E_LMPC, delimiter=',', fmt='%f')99

    if LMPC.Solver == "CVX":
        G_LMPC_sparse = spmatrix(G_LMPC[np.nonzero(G_LMPC)], np.nonzero(G_LMPC)[0].astype(int),
                                 np.nonzero(G_LMPC)[1].astype(int), G_LMPC.shape)
        E_LMPC_sparse = spmatrix(E_LMPC[np.nonzero(E_LMPC)], np.nonzero(E_LMPC)[0].astype(int),
                                 np.nonzero(E_LMPC)[1].astype(int), E_LMPC.shape)
        G_LMPC_return = G_LMPC_sparse
        E_LMPC_return = E_LMPC_sparse
    else:
        G_LMPC_return = G_LMPC
        E_LMPC_return = E_LMPC

    return G_LMPC_return, E_LMPC_return


# 等式约束
def _LMPC_BuildMatEqConst(LMPC, A, B, C, N, n, d):
    # Buil matrices for optimization (Convention from Chapter 15.2 Borrelli, Bemporad and Morari MPC book)
    # We are going to build our optimization vector z \in \mathbb{R}^((N+1) \dot n \dot N \dot d), note that this vector
    # stucks the predicted trajectory x_{k|t} \forall k = t, \ldots, t+N+1 over the horizon and
    # the predicte input u_{k|t} \forall k = t, \ldots, t+N over the horizon
    Gx = np.eye(n * (N + 1))
    Gu = np.zeros((n * (N + 1), d * (N)))

    E = np.zeros((n * (N + 1), n))
    E[np.arange(n)] = np.eye(n)

    L = np.zeros((n * (N + 1) + n + 1, 1))  # n+1 for the terminal constraint
    L[-1] = 1  # Summmation of lamba must add up to 1

    for i in range(0, N):
        ind1 = n + i * n + np.arange(n)
        ind2x = i * n + np.arange(n)
        ind2u = i * d + np.arange(d)

        Gx[np.ix_(ind1, ind2x)] = -A[i]
        Gu[np.ix_(ind1, ind2u)] = -B[i]
        L[ind1, :] = C[i]

    G = np.hstack((Gx, Gu))

    if LMPC.Solver == "CVX":
        L_sparse = spmatrix(L[np.nonzero(L)], np.nonzero(L)[0].astype(int), np.nonzero(L)[1].astype(int), L.shape)
        L_return = L_sparse
    else:
        L_return = L

    return L_return, G, E


def _LMPC_GetPred(Solution, n, d, N, np):
    xPred = np.squeeze(np.transpose(np.reshape((Solution[np.arange(n * (N + 1))]), (N + 1, n))))
    uPred = np.squeeze(np.transpose(np.reshape((Solution[n * (N + 1) + np.arange(d * N)]), (N, d))))
    lambd = Solution[(n * (N + 1) + d * N):(Solution.shape[0] - n - 2 * N)]
    slack = Solution[Solution.shape[0] - n - 2 * N:Solution.shape[0] - 2 * N]
    laneSlack = Solution[Solution.shape[0] - 2 * N:]

    # print np.sum(np.abs(laneSlack))
    # if np.sum(np.abs(laneSlack)) > 0.5:
    #     pdb.set_trace()

    return xPred, uPred, lambd, slack


# ======================================================================================================================
# ======================================================================================================================
# ========================= Internal functions for Local Regression and Linearization ==================================
# ======================================================================================================================
# ======================================================================================================================
def  _LMPC_EstimateABC(ControllerLMPC, sortedLapTime):
    LinPoints = ControllerLMPC.LinPoints
    LinInput = ControllerLMPC.LinInput
    N = ControllerLMPC.N
    n = ControllerLMPC.n
    d = ControllerLMPC.d
    TimeSS = ControllerLMPC.TimeSS
    LapCounter = ControllerLMPC.LapCounter
    PointAndTangent = ControllerLMPC.map.PointAndTangent
    dt = ControllerLMPC.dt
    it = ControllerLMPC.it
    SS = ControllerLMPC.SS
    uSS = ControllerLMPC.uSS
    tSSparaTime = ControllerLMPC.tSSparaTime


    ParallelComputation = 0
    Atv = [];
    Btv = [];
    Ctv = [];
    indexUsed_list = []

    usedIt = sortedLapTime[
             0:ControllerLMPC.itUsedSysID]  # range(ControllerLMPC.it-ControllerLMPC.itUsedSysID, ControllerLMPC.it)是挑选的圈数，默认的itusedsysid=1
    MaxNumPoint = 40  # Need to reason on how these points are selected

    # 将预测时域内每个时刻的Ai,Bi,Ci,都添加到Atv，Btv，Ctv中
    for i in range(0, N):
        # 最后传入tSSparaTime参数，用来删除每一小段的最后一个点，避免最后一个点没有下一个点
    #    print 'n=%d'%(i)


        Ai, Bi, Ci, indexSelected = RegressionAndLinearization(LinPoints, LinInput, usedIt, SS, uSS, TimeSS,
                                                               MaxNumPoint, qp, n, d, matrix, PointAndTangent, dt, i,tSSparaTime)
        Atv.append(Ai)
        Btv.append(Bi)
        Ctv.append(Ci)
        indexUsed_list.append(indexSelected)

    return Atv, Btv, Ctv, indexUsed_list


def RegressionAndLinearization(LinPoints, LinInput, usedIt, SS, uSS, LapCounter, MaxNumPoint, qp, n, d, matrix,
                               PointAndTangent, dt, i,tSSparaTime):
    # i是滚动时域N中遍历的第i个点

    x0 = LinPoints[i, :]

    Ai = np.zeros((n, n))
    Bi = np.zeros((n, d))
    Ci = np.zeros((n, 1))

    # Compute Index to use
    h = 5
    lamb = 0.0
    stateFeatures = [0, 1, 2]
    ConsiderInput = 1

    if ConsiderInput == 1:
        scaling = np.array([[0.1, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 1.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 1.0]])
        xLin = np.hstack((LinPoints[i, stateFeatures], LinInput[i, :]))
    else:
        scaling = np.array([[1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0]])
        xLin = LinPoints[i, stateFeatures]

    indexSelected = []
    K = []
    for ii in usedIt:
        indexSelected_i, K_i = ComputeIndex(h, SS, uSS, LapCounter, ii, xLin, stateFeatures, scaling, MaxNumPoint,
                                            ConsiderInput, tSSparaTime[:,ii])

        indexSelected.append(indexSelected_i)
        K.append(K_i)

    # =========================
    # ====== Identify vx ======
    inputFeatures = [1]
    Q_vx, M_vx = Compute_Q_M(SS, uSS, indexSelected, stateFeatures, inputFeatures, usedIt, np, matrix, lamb, K)

    yIndex = 0
    b = Compute_b(SS, yIndex, usedIt, matrix, M_vx, indexSelected, K, np)
    Ai[yIndex, stateFeatures], Bi[yIndex, inputFeatures], Ci[yIndex] = LMPC_LocLinReg(Q_vx, b, stateFeatures,
                                                                                        inputFeatures, qp)

    # =======================================
    # ====== Identify Lateral Dynamics ======
    inputFeatures = [0]
    Q_lat, M_lat = Compute_Q_M(SS, uSS, indexSelected, stateFeatures, inputFeatures, usedIt, np, matrix, lamb, K)

    yIndex = 1  # vy
    b = Compute_b(SS, yIndex, usedIt, matrix, M_lat, indexSelected, K, np)
    Ai[yIndex, stateFeatures], Bi[yIndex, inputFeatures], Ci[yIndex] = LMPC_LocLinReg(Q_lat, b, stateFeatures,
                                                                                      inputFeatures, qp)

    yIndex = 2  # wz
    b = Compute_b(SS, yIndex, usedIt, matrix, M_lat, indexSelected, K, np)
    Ai[yIndex, stateFeatures], Bi[yIndex, inputFeatures], Ci[yIndex] = LMPC_LocLinReg(Q_lat, b, stateFeatures,
                                                                                      inputFeatures, qp)
    # 在速度变量x0[0,1,2]的情况下，位置变量x0[3,4,5]的变化
    # ===========================
    # ===== Linearization =======
    vx = x0[0];
    vy = x0[1]
    wz = x0[2];
    epsi = x0[3]
    s = x0[4];
    ey = x0[5]

    if s < 0:
        print "s is negative, here the state: \n", LinPoints

    startTimer = datetime.datetime.now()  # Start timer for LMPC iteration
    cur = Curvature(s, PointAndTangent)
    cur = Curvature(s, PointAndTangent)
    den = 1 - cur * ey
    # ===========================
    # ===== Linearize epsi ======
    # epsi_{k+1} = epsi + dt * ( wz - (vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey) * cur )
    depsi_vx = -dt * np.cos(epsi) / den * cur
    depsi_vy = dt * np.sin(epsi) / den * cur
    depsi_wz = dt
    depsi_epsi = 1 - dt * (-vx * np.sin(epsi) - vy * np.cos(epsi)) / den * cur
    depsi_s = 0  # Because cur = constant
    depsi_ey = dt * (vx * np.cos(epsi) - vy * np.sin(epsi)) / (den ** 2) * cur * (-cur)

    Ai[3, :] = [depsi_vx, depsi_vy, depsi_wz, depsi_epsi, depsi_s, depsi_ey]
    Ci[3] = epsi + dt * (wz - (vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey) * cur) - np.dot(Ai[3, :], x0)
    # ===========================
    # ===== Linearize s =========
    # s_{k+1} = s    + dt * ( (vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey) )
    ds_vx = dt * (np.cos(epsi) / den)
    ds_vy = -dt * (np.sin(epsi) / den)
    ds_wz = 0
    ds_epsi = dt * (-vx * np.sin(epsi) - vy * np.cos(epsi)) / den
    ds_s = 1  # + Ts * (Vx * cos(epsi) - Vy * sin(epsi)) / (1 - ey * rho) ^ 2 * (-ey * drho);
    ds_ey = -dt * (vx * np.cos(epsi) - vy * np.sin(epsi)) / (den * 2) * (-cur)

    Ai[4, :] = [ds_vx, ds_vy, ds_wz, ds_epsi, ds_s, ds_ey]
    Ci[4] = s + dt * ((vx * np.cos(epsi) - vy * np.sin(epsi)) / (1 - cur * ey)) - np.dot(Ai[4, :], x0)

    # ===========================
    # ===== Linearize ey ========
    # ey_{k+1} = ey + dt * (vx * np.sin(epsi) + vy * np.cos(epsi))
    dey_vx = dt * np.sin(epsi)
    dey_vy = dt * np.cos(epsi)
    dey_wz = 0
    dey_epsi = dt * (vx * np.cos(epsi) - vy * np.sin(epsi))
    dey_s = 0
    dey_ey = 1

    Ai[5, :] = [dey_vx, dey_vy, dey_wz, dey_epsi, dey_s, dey_ey]
    Ci[5] = ey + dt * (vx * np.sin(epsi) + vy * np.cos(epsi)) - np.dot(Ai[5, :], x0)

    endTimer = datetime.datetime.now();
    deltaTimer_tv = endTimer - startTimer

    return Ai, Bi, Ci, indexSelected


def Compute_Q_M(SS, uSS, indexSelected, stateFeatures, inputFeatures, usedIt, np, matrix, lamb, K):
    Counter = 0
    it = 1
    X0 = np.empty((0, len(stateFeatures) + len(inputFeatures)))
    Ktot = np.empty((0))

    for it in usedIt:
        X0 = np.append(X0, np.hstack((np.squeeze(SS[np.ix_(indexSelected[Counter], stateFeatures, [it])]),
                                      np.squeeze(uSS[np.ix_(indexSelected[Counter], inputFeatures, [it])], axis=2))),
                       axis=0)
        Ktot = np.append(Ktot, K[Counter])
        Counter = Counter + 1
    # 横向拼接
    M = np.hstack((X0, np.ones((X0.shape[0], 1))))
    Q0 = np.dot(np.dot(M.T, np.diag(Ktot)), M)
    Q = matrix(Q0 + lamb * np.eye(Q0.shape[0]))

    return Q, M


# b=-M.T*K*Y
def Compute_b(SS, yIndex, usedIt, matrix, M, indexSelected, K, np):
    Counter = 0
    y = np.empty((0))
    Ktot = np.empty((0))

    for it in usedIt:
        y = np.append(y, np.squeeze(SS[np.ix_(indexSelected[Counter] + 1, [yIndex], [it])]))
        Ktot = np.append(Ktot, K[Counter])
        Counter = Counter + 1

    b = matrix(-np.dot(np.dot(M.T, np.diag(Ktot)), y))

    return b


def LMPC_LocLinReg(Q, b, stateFeatures, inputFeatures, qp):
    import numpy as np
    from numpy import linalg as la
    import datetime

    # K = np.ones(len(index))
    # 求解令(1/2)（X.T*Q*X+b*X）最小的X
    startTimer = datetime.datetime.now()  # Start timer for LMPC iteration
    res_cons = qp(Q, b)  # This is ordered as [A B C]



    endTimer = datetime.datetime.now();
    deltaTimer_tv = endTimer - startTimer

    # print "Non removable time: ", deltaTimer_tv.total_seconds()
    Result = np.squeeze(np.array(res_cons['x']))
    A = Result[0:len(stateFeatures)]
    B = Result[len(stateFeatures):(len(stateFeatures) + len(inputFeatures))]
    C = Result[-1]

    return A, B, C


def ComputeIndex(h, SS, uSS, LapCounter, it, x0, stateFeatures, scaling, MaxNumPoint, ConsiderInput,tSSparaTime):
    import numpy as np
    from numpy import linalg as la
    import datetime

    startTimer = datetime.datetime.now()  # Start timer for LMPC iteration

    # What to learn a model such that: x_{k+1} = A x_k  + B u_k + C
    oneVec = np.ones((SS[0:LapCounter[it], :, it].shape[0] - 1, 1))
    # 一列x0,横着是x0的六个状态向量
    x0Vec = (np.dot(np.array([x0]).T, oneVec.T)).T
    # 考虑输入就把输入变量和状态变量放一块
    if ConsiderInput == 1:
        DataMatrix = np.hstack((SS[0:LapCounter[it] - 1, stateFeatures, it], uSS[0:LapCounter[it] - 1, :, it]))
        # 没有输入就只要状态变量
    else:
        DataMatrix = SS[0:LapCounter[it] - 1, stateFeatures, it]

    # np.savetxt('A.csv', SS[0:TimeSS[it]-1, stateFeatures, it], delimiter=',', fmt='%f')
    # np.savetxt('B.csv', SS[0:TimeSS[it], :, it][0:-1, stateFeatures], delimiter=',', fmt='%f')
    # np.savetxt('SS.csv', SS[0:TimeSS[it], :, it], delimiter=',', fmt='%f')
    # 挑选出和当前状态最近的点
    diff = np.dot((DataMatrix - x0Vec), scaling)
    norm = la.norm(diff, 1, axis=1)
    indexdell=np.where(tSSparaTime.astype(int) - norm.shape>=0)[0][0]
    norm[tSSparaTime.astype(int)[0:indexdell]-1]=1000 #将每小段最后一个点和当前点的距离设置为特别大，将其排除
    indexTot = np.squeeze(np.where(norm < h))
  #  if any(indexTot)==False:

    if indexTot.shape<10:
        print 'indexTot', indexTot.shape



    # print indexTot.shape, np.argmin(norm), norm, x0
    # 如果和当前点很近的点过多，只选取maxnumpoint个点
    if (indexTot.shape >= MaxNumPoint):
        index = np.argsort(norm)[0:MaxNumPoint]
        # MinNorm = np.argmin(norm)
        # if MinNorm+MaxNumPoint >= indexTot.shape[0]:
        #     index = indexTot[indexTot.shape[0]-MaxNumPoint:indexTot.shape[0]]
        # else:
        #     index = indexTot[MinNorm:MinNorm+MaxNumPoint]
    else:
        index = indexTot


    K = (1 - (norm[index] / h) ** 2) * 3 / 4
    # K = np.ones(len(index))

    return index, K