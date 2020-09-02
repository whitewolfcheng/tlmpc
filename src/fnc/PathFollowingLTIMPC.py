# coding=UTF-8
from scipy import linalg
import numpy as np
from cvxopt.solvers import qp
from cvxopt import spmatrix, matrix, solvers
import datetime, pdb
#使用cvxopt求解器求解凸优化问题，凸优化问题描述：
# min (X^T*M*X+q^T*X)/2
#s.t.Fx≤b，Gx=E
#cvxopt求解器：cvxopt.solvers.qp(P,q,G,h,A,b)
solvers.options['show_progress'] = False

class PathFollowingLTI_MPC:
    """Create the Path Following LMPC controller with LTI model
    Attributes:
        solve: given x0 computes the control action
    """
    def __init__(self, A, B, Q, R, N, vt, inputConstr):
        """Initialization
        A, B: Liner Time Invariant (LTI) system dynamics
        Q, R: weights to build the cost function h(x,u) = ||x||_Q + ||u||_R
        N: horizon length
        vt: target velocity
        """
        self.A = A         # n x n matrix
        self.B = B         # n x d matrix
        self.n = A.shape[0]   # number of states
        self.d = B.shape[1]   # number of inputs
        self.N = N            # horizon length
        self.Q = Q            #代价函数系数 h(x,u) = ||x||_Q + ||u||_R
        self.R = R
        self.vt = vt        # target velocity
        startTimer = datetime.datetime.now()
        endTimer = datetime.datetime.now(); deltaTimer = endTimer - startTimer
        self.solverTime = deltaTimer
        self.linearizationTime = deltaTimer
        self.inputConstr = inputConstr

        self.M, self.q = _buildMatCost(self)
        self.F, self.b = _buildMatIneqConst(self)
        self.G, self.E = _buildMatEqConst(self)

    def solve(self, x0):
        M = self.M
        F = self.F
        G = self.G
        E = self.E
        q = self.q
        b = self.b
        n, d = self.n, self.d
        N = self.N

        startTimer = datetime.datetime.now()
        sol = qp(M, matrix(q), F, matrix(b), G, E * matrix(x0)) #代价函数系数min (X^T*M*X+q^T*X)/2、不等式约束系数Fx≤b、等式约束系数Gx=E

        endTimer = datetime.datetime.now(); deltaTimer = endTimer - startTimer
        self.solverTime = deltaTimer
        if sol['status'] == 'optimal':
            self.feasible = 1
        else:
            self.feasible = 0

        self.xPred = np.squeeze(np.transpose(np.reshape((np.squeeze(sol['x'])[np.arange(n * (N + 1))]), (N + 1, n)))).T
        #transpose是转置，两次转置再穿插reshape和删除一维条目squeeze，将(N+1)*n个[x]元素变成n组向量，每个向量代表一种状态变量，向量中包含它的N+1个数据
        self.uPred = np.squeeze(np.transpose(np.reshape((np.squeeze(sol['x'])[n * (N + 1) + np.arange(d * N)]), (N, d)))).T



# ======================================================================================================================
# ======================================================================================================================
# =============================== Internal functions for MPC reformulation to QP =======================================
# ======================================================================================================================
# ======================================================================================================================
#求解等式约束系数Gx=E
def _buildMatEqConst(Controller):
    # Buil matrices for optimization (Convention from Chapter 15.2 Borrelli, Bemporad and Morari MPC book)
    # We are going to build our optimization vector z \in \mathbb{R}^((N+1) \dot n \dot N \dot d), note that this vector
    # stucks the predicted trajectory x_{k|t} \forall k = t, \ldots, t+N+1 over the horizon and
    # the predicte input u_{k|t} \forall k = t, \ldots, t+N over the horizon
    n = Controller.n #状态变量数量
    N = Controller.N #滚动长度
    d = Controller.d #输入变量数量
    A,B = Controller.A, Controller.B

    Gx = np.eye(n * (N + 1))
    Gu = np.zeros((n * (N + 1), d * (N)))

    E = np.zeros((n * (N + 1), n))
    E[np.arange(n)] = np.eye(n)   #上面是一个n*n的单位矩阵，下面是零矩阵

    for i in range(0, N):#当前时刻状态变量和控制变量的系数，乘以A、B使得当前时刻变成下一时刻的估计值，加上负号使得下一时刻的估计值与右边实际的下一时刻值相减，最后等于0
        ind1 = n + i * n + np.arange(n) #从第二组开始，每次选够x的长度个坐标
        ind2x = i * n + np.arange(n)
        ind2u = i * d + np.arange(d)
       #np.ix 花式索引，可选取在每个ind1行中的ind2x列，最后得出size.ind1行，size.ind2列的一个矩阵
        Gx[np.ix_(ind1, ind2x)] = -A
        Gu[np.ix_(ind1, ind2u)] = -B

    G = np.hstack((Gx, Gu))

    G_sparse = spmatrix(G[np.nonzero(G)], np.nonzero(G)[0].astype(int), np.nonzero(G)[1].astype(int), G.shape)
    E_sparse = spmatrix(E[np.nonzero(E)], np.nonzero(E)[0].astype(int), np.nonzero(E)[1].astype(int), E.shape)

    return G_sparse, E_sparse
#求解不等式约束系数Fx≤b（状态变量和输入变量的约束）
def _buildMatIneqConst(Controller):
    N = Controller.N
    n = Controller.n
    inputConstr = Controller.inputConstr
#状态约束的矩阵
    # Buil the matrices for the state constraint in each region. In the region i we want Fx[i]x <= bx[b]
    Fx = np.array([[1., 0., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 0., 1.],
                   [0., 0., 0., 0., 0., -1.]])

    bx = np.array([[10],  # vx max
                   [2.],  # max ey
                   [2.]]), # max ey
#输入约束的矩阵
    # Buil the matrices for the input constraint in each region. In the region i we want Fx[i]x <= bx[b]
    Fu = np.array([[1., 0.],
                   [-1., 0.],
                   [0., 1.],
                   [0., -1.]])

    bu = np.array([[inputConstr[0,0]],  # -Min Steering
                   [inputConstr[0,1]],  # Max Steering
                   [inputConstr[1,0]],  # -Min Acceleration
                   [inputConstr[1,1]]])  # Max Acceleration
#整合输入约束和状态约束矩阵
    #状态约束矩阵
    # Now stuck the constraint matrices to express them in the form Fz<=b. Note that z collects states and inputs
    # Let's start by computing the submatrix of F relates with the state
    rep_a = [Fx] * (N)
    Mat = linalg.block_diag(*rep_a)
    NoTerminalConstr = np.zeros((np.shape(Mat)[0], n))  # No need to constraint also the terminal point 终端约束，因为不需要终端约束所以设置为0矩阵
    Fxtot = np.hstack((Mat, NoTerminalConstr))
    bxtot = np.tile(np.squeeze(bx), N)   #删除bx数组形状中的单维数组条目然后沿x方向复制N遍（其实应该是列向量，但求解器中的输入是行向量）
    #输入约束矩阵
    # Let's start by computing the submatrix of F relates with the input
    rep_b = [Fu] * (N)
    Futot = linalg.block_diag(*rep_b)
    butot = np.tile(np.squeeze(bu), N)
    #整合输入和状态约束矩阵
    # Let's stack all together
    rFxtot, cFxtot = np.shape(Fxtot)
    rFutot, cFutot = np.shape(Futot)
    Dummy1 = np.hstack((Fxtot, np.zeros((rFxtot, cFutot))))
    Dummy2 = np.hstack((np.zeros((rFutot, cFxtot)), Futot))
    F = np.vstack((Dummy1, Dummy2))
    b = np.hstack((bxtot, butot))

    F_sparse = spmatrix(F[np.nonzero(F)], np.nonzero(F)[0].astype(int), np.nonzero(F)[1].astype(int), F.shape)
    return F_sparse, b
#求解代价函数系数min (X^T*M*X+q*X)/2
def _buildMatCost(Controller):
    N = Controller.N    #滚动长度
    Q = Controller.Q    #代价函数系数 h(x,u) = ||x||_Q + ||u||_R
    R = Controller.R
    vt = Controller.vt  #目标速度
    P = Controller.Q
    b = [Q] * (N)
    #linalg.block_diag:创建广义主对角线非零大矩阵，输入参数是矩阵，将矩阵作为广义主对角线元素
    Mx = linalg.block_diag(*b)

    c = [R] * (N)
    Mu = linalg.block_diag(*c)

    M0 = linalg.block_diag(Mx, P, Mu)
    xtrack = np.array([vt, 0, 0, 0, 0, 0])
    # np.tile将xtrack沿着x轴复制N+1遍，np.append将两个数组横向合并
    q = - 2 * np.dot(np.append(np.tile(xtrack, N + 1), np.zeros(R.shape[0] * N)), M0)  #代表给定速度的量，为负数，实际速度减去给定速度作为速度误差
    M = 2 * M0  # Need to multiply by two because CVX considers 1/2 in front of quadratic cost
    # np.nonzero函数提取矩阵中非零元素的坐标,spmatrix创建稀疏矩阵,经实验后发现spmatrix操作后和原矩阵一样，只是里面的元素变成了1e00这种格式，应该是变成稀疏矩阵的格式方便计算提高精度
    M_sparse = spmatrix(M[np.nonzero(M)], np.nonzero(M)[0].astype(int), np.nonzero(M)[1].astype(int), M.shape)
    return M_sparse, q