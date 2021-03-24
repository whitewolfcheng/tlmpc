#coding=utf-8
import numpy as np
from numpy import linalg as la
from Utilities import Curvature
from numpy import hstack
from scipy import interpolate #插值

class ControllerTLMPC():
    """Create the LMPC
    Attributes:
        solve: given x0 computes the control action
        addTrajectory: given a ClosedLoopData object adds the trajectory to SS, Qfun, uSS and updates the iteration index
        addPoint: this function allows to add the closed loop data at iteration j to the SS of iteration (j-1)
        update: this function can be used to set SS, Qfun, uSS and the iteration index.
    """

    def __init__(self, numSS_Points,numSR_Points, N,  dt,  map, Laps,TLaps, TimeLMPC ):
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
        self.numSR_points=numSR_Points  #similar road 即记忆中和当前道路相似的轨迹
        self.N = N
        self.n = 6
        self.d = 2
        self.dt = dt
        self.map = map
        self.LapTime = 0
        self.TLaps = TLaps

        self.OldInput = np.zeros((1,2))
#初始化
        # Initialize the following quantities to avoid dynamic allocation
        NumPoints = int(TimeLMPC / dt)*2 + 1
        self.LapCounter = 10000 * np.ones(Laps).astype(int)        # Time at which each j-th iteration is completed
        self.TimeSS     = 10000 * np.ones(Laps).astype(int)        # Time at which each j-th iteration is completed
        self.transSSTime  = 0*np.ones(NumPoints).astype(int)  #迁移后的SS集中的点个数(分段统计)
        self.tSSlapsTime= 0*np.ones(TLaps).astype(int)  #迁移后的SS集中每一圈的点个数
        self.tSSparaTime= 0*np.ones((NumPoints,TLaps)) #迁移后，每一圈里各片段的点个数
        self.SS         = 10000 * np.ones((NumPoints, 6, Laps))    # Sampled Safe SS
        self.transSS    = 10000 * np.ones((NumPoints,6))
        self.transuSS  = 10000 * np.ones((NumPoints, 2))
        self.BestbeforePoint = 10000 * np.ones((NumPoints, 6))
        self.chazhi_PointSelect = 10000 * np.ones((NumPoints, 6))
        self.uSS        = 10000 * np.ones((NumPoints, 2, Laps))    # Input associated with the points in SS
        self.bcur       = 100000 * np.ones((NumPoints, Laps))       # bcur: 数据集中各点处的曲率
        self.arclen = 0    * np.ones((NumPoints, Laps))            # arclen:该点与前一点之间的弧长
        self.arclens = 0 * np.ones((NumPoints, Laps))               # arclens:起点到该点间的弧长
        self.SS_glob    = 10000 * np.ones((NumPoints, 6, Laps))    # SS in global (X-Y) used for plotting
        self.overflag=0
        self.overlap=0
        self.passway=0
        self.TranslapTime=0

        self.zVector = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) #总路程19.22

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
        SS = self.SS;
      #  SS_glob=self.SS_glob
        bcur = self.bcur
        arclen=self.arclen
        arclens=self.arclens
        uSS = self.uSS;
        TimeSS = self.TimeSS
        LapCounter=self.LapCounter
        OldInput = self.OldInput
        N = self.N;
        dt = self.dt
        it = self.it
        numSS_Points = self.numSS_Points
        numSR_Points=self.numSR_points
        LinPoints = self.LinPoints
        LinInput = self.LinInput
        map = self.map
        PointAndTangent = map.PointAndTangent

        # Select Points from SS

        #    self.LinPoints[4, -1] = self.LinPoints[4, -1] - map.TrackLength
        #判断相似路径
        SR_PointSelectedTot = np.empty((n, 0))

        Succ_BSS_PointSelectedTot = np.empty((n, 0))
        Succ_BuSS_PointSelectedTot = np.empty((d, 0))
        Qfun_SelectedTot = np.empty((0))
        reallyroadarclen=np.empty(numSR_Points)
        maxsspoint =max(TimeSS)
        Mindiff=100.0 #先整个大的
        jdelerror=[] #由于删除长度不够的列，导致的jdel与实际MINJJ的误差
        MinJarray = range(0, 42)
        for ii in range(0,maxsspoint):
            #剔除长度不够的列
            jdell=2
            for jdel in range(2,42):#jdel是真正的圈次数
                if jdell >= bcur.shape[1] :
                    break
                while jdel in jdelerror:
                    jdel+=1          #有可能jdel+1也在jdelerror里面，因此需要while
                if ii+numSR_Points>TimeSS[jdel]-1 or ii==LapCounter[jdel] :      #注意最后一个点没办法计算曲率，因此实际的点个数应少1
                    bcur=np.delete(bcur,jdell,axis=1)
                    arclens=np.delete(arclens,jdell,axis=1)
                    jdell=jdell-1
                    if jdel not in  jdelerror:
                        jdelerror.append(jdel)
                        MinJarray.remove(jdel)

                jdell=jdell+1

            #记忆中的numSR_Points个点的曲率
            curb=bcur[ii:ii+numSR_Points,2:]

            #当前路径中心线的numSR_Points个点的对应s值及曲率
            nowarclen=arclens[ii:ii+numSR_Points,2:]  #s值

            if curb.shape[1]==0: #如果曲率是空集，则跳出循环
                continue

            nowxs=self.zVector[4]  #修正s值
       #     diffarclens=nowxs-nowarclen[9,:]
        #    nowarclens=nowarclen+diffarclens
            #如果在最开始，zvector[4]很小，使得nowarclens中出现负数，则片段路径的需要以第一个点为基准点
         #   if np.any(nowarclens<0):
            if nowarclen[-1,-1]==0:  #最后一行最后一列会出现长度是0的情况
                break
            diffarclens = nowxs - nowarclen[0, :]
            nowarclens = nowarclen + diffarclens
            #利用弧长求曲率
            nowcur = np.empty((nowarclen.shape[0], nowarclen.shape[1]))
            for i in range(0,nowarclen.shape[0]):
                for j in range(0, nowarclen.shape[1]):
                    if nowarclens[i,j]<0:
                       # print nowarclens[i,j]+map.TrackLength
                        nowcur[i, j] = Curvature(nowarclens[i, j]+map.TrackLength, PointAndTangent)
                    else:
                        nowcur[i,j] = Curvature(nowarclens[i,j], PointAndTangent)
            #现实路径与过往轨迹的曲率差的范数
            diffcur=nowcur-curb
            diffcurnorm = la.norm(diffcur, 2, axis=0)
            MinJ = np.argmin(diffcurnorm)
            subMinarclens=nowarclens[:,MinJ]  #待定的最小路程
            subnowcur=nowcur[:,MinJ]
            #将从被删除后的矩阵中选取的MinJ变成完整矩阵中的MinJ
            if Mindiff>diffcurnorm[MinJ]:
                MinJJ=MinJ
                MinII=ii
                Mindiff=diffcurnorm[MinJ]

                Minarclens=subMinarclens
                MinJJ=MinJarray[MinJJ+2] #原矩阵的索引值
                minnowcur=subnowcur

        # 虚拟路径跟踪轨迹
        similarroadpoint=SS[MinII:MinII+numSR_Points,:, MinJJ]
        self.TSR_PointSelectedTot = similarroadpoint

        # 现实道路中心线
        reallyroadarclen = Minarclens.reshape((numSR_Points, 1))
        self.reallyroadPoint = hstack((hstack((np.zeros((reallyroadarclen.shape[0], 4)), reallyroadarclen)),
                                       np.zeros((reallyroadarclen.shape[0], 1))))  # 实际道路为只有s量，其余各量为0的点
        #对虚拟路径跟踪轨迹进行插值，先求样条参数
        similarroadpoint_s = SS[MinII:MinII + numSR_Points, 4, MinJJ] #以s值为参考进行插值
        similarroadpoint_Vx = SS[MinII:MinII + numSR_Points, 0, MinJJ]
        similarroadpoint_Vy = SS[MinII:MinII + numSR_Points, 1, MinJJ]
        similarroadpoint_psi = SS[MinII:MinII + numSR_Points, 3, MinJJ] #先找到虚拟路径中心线的各参数
        similarroadpoint_wz = SS[MinII:MinII + numSR_Points, 2, MinJJ]

        similarroadpoint_arclens = Minarclens.reshape((1, numSR_Points))  #arclens其实就是当前坐标系下的s值
        similarroadpoint_ey = SS[MinII:MinII + numSR_Points, 5, MinJJ]

        tck_Vx = interpolate.interp1d(similarroadpoint_s, similarroadpoint_Vx)
        tck_Vy = interpolate.interp1d(similarroadpoint_s, similarroadpoint_Vy)
        tck_wz = interpolate.interp1d(similarroadpoint_s, similarroadpoint_wz)
        tck_psi = interpolate.interp1d(similarroadpoint_s, similarroadpoint_psi,'cubic')  #虚拟路径跟踪轨迹的样条插值参数
        tck_arclens = interpolate.interp1d(similarroadpoint_s, similarroadpoint_arclens,3 )
        tck_y= interpolate.interp1d(similarroadpoint_s, similarroadpoint_ey,3)


           #开始经验迁移，找到最优轨迹
        Sorigin=SS[MinII,4,MinJJ]  #相似道路中心线起点的s值
        Send = SS[MinII+numSR_Points-1, 4, MinJJ] #相似道路中心线终点的s值
        MinNumPoints=10000 #先整个大的初始值
        MinNumerror=10000
        #上一段最优轨迹的最后一个点
        if self.transSSTime[self.TransTime-1]==0:
            BestzVector=self.zVector
            BestzVectors =self.zVector
        else :
            BestzVector = self.transSS[self.transSSTime[self.TransTime - 1] - 1, :]
            BestzVectors= self.transSS[self.transSSTime[self.TransTime-1]-11:self.transSSTime[self.TransTime-1], :]

        for J in range(SS.shape[2]-1,1,-1): #为了快点找到最优轨迹，倒着遍历
            totalpointS=SS[ : ,4, J ] #所有点的s值
            TrajectoriesIndex=np.argwhere((totalpointS<=Send) & (totalpointS>=Sorigin)) #s值在相似道路中心线范围内的轨迹点索引值
            selectedTrajecties = SS[TrajectoriesIndex[:,0], :, J]
            difflastpoint = selectedTrajecties-BestzVector
            normselectedtrajecties = la.norm(difflastpoint, 1, axis=1)
            MinNormselectedtrajecties = np.argmin(normselectedtrajecties)
            if MinNormselectedtrajecties-10<0:
                similarLastpoint = selectedTrajecties[MinNormselectedtrajecties, :]
            else:
                similarLastpoint=selectedTrajecties[MinNormselectedtrajecties-10:MinNormselectedtrajecties+1,:]
        #    if selectedTrajecties.shape[0]>numSR_Points: #待选轨迹不能比沿着路径中心线跑还要慢
         #       continue
        #    print la.norm(self.zVector[0:3]-selectedTrajecties[0,0:3],1)

            if selectedTrajecties.shape[0]+la.norm(similarLastpoint-BestzVectors,1).sum()*10>MinNumerror: #待选轨迹比最优轨迹还慢的也扔了吧,同时保证速度变化慢
                continue
            if selectedTrajecties.shape[0]<numSS_Points/2+1: #待选轨迹如果太快，会导致点数过少无法建立ss集，加1是为了划分当前SS集与下一步的SS集
                continue
            #X1表示待迁移的轨迹点，X2表示虚拟道路中心线上的点
            TransX1=selectedTrajecties[:,:]
            TransX2_s=TransX1[:,4]

            #找到X2
         #   TransX2_Vx = tck_Vx(TransX2_s)
         #   TransX2_Vy = tck_Vy(TransX2_s)
          #  TransX2_WZ = tck_wz(TransX2_s)
            TransX2_epsi = tck_psi(TransX2_s)
            TransX2_arclens=tck_arclens(TransX2_s)
            TransX2_ey = tck_y(TransX2_s)
            #迁移X1
            TransX1_ey=(TransX1[:,5]-TransX2_ey)*np.cos(TransX2_epsi)
            if max(TransX1_ey)>map.halfWidth:  #删除大于路宽约束的点
                continue
            TransX1_s=TransX2_arclens+(TransX1[:,5]-TransX2_ey)*np.sin(TransX2_epsi) #先是按原来轨迹计算的s值你
            TransX1_epsi=TransX1[:,3]-TransX2_epsi
            # X1处输入变量
            TransU1 = uSS[TrajectoriesIndex[:, 0], :, J]
       #     TransX1_WZ=TransX1[:,2]
        #    TransX1_Vy=TransX1[:,1]*np.cos(TransX2_epsi)-TransX1[:, 0] * np.sin(TransX2_epsi)
        #    TransX1_Vx = TransX1[:, 0] * np.cos(TransX2_epsi)+TransX1[:, 1] * np.sin(TransX2_epsi)

            MinNumPoints=selectedTrajecties.shape[0]
            MinNumerror=MinNumPoints+la.norm(similarLastpoint-BestzVectors,1).sum()*10
            BestSS_PointSelect = np.empty((n, MinNumPoints))
            BestSS_PointSelect[0:3, :] = TransX1[:,0:3].T
            BestSS_PointSelect[3, :] = TransX1_epsi
            BestSS_PointSelect[4, :] = TransX1_s
            BestSS_PointSelect[5, :] = TransX1_ey
            #检验插值
            chazhi_PointSelect = np.empty((n, MinNumPoints))
            chazhi_PointSelect[3, :] = TransX2_epsi
            chazhi_PointSelect[4, :] = TransX2_arclens
            chazhi_PointSelect[5, :] = TransX2_ey
            BestbeforePoint=TransX1 #迁移前的轨迹

        TransTime=self.TransTime #仿真的第几段路段，在Virtualsimulater程序中每次加1
        transSSTime=self.transSSTime[TransTime-1] #该路段内最优轨迹总共有几个点，self.transSSTime是统计所有transSSTime的集合
    #    if transSSTime+MinNumPoints >= self.transSS.shape[0] :
    #        self.transSS[transSSTime:, :] = BestSS_PointSelect[:, :self.transSS.shape[0]-transSSTime].T   #防止超过最大仿真时间
     #       self.BestbeforePoint[transSSTime: , :] = BestbeforePoint[:, :self.transSS.shape[0]-transSSTime]
     #       self.chazhi_PointSelect[transSSTime:, :] = chazhi_PointSelect[:, :self.transSS.shape[0]-transSSTime].T
     #       self.overflag=1 #超出最大仿真时间了

        self.transSS[transSSTime:transSSTime+MinNumPoints,:] =  BestSS_PointSelect.T  #将迁移后的最优轨迹加入新的ss集

        self.transuSS[transSSTime:transSSTime+MinNumPoints,:]= TransU1
        self.BestbeforePoint[transSSTime:transSSTime + MinNumPoints, :] = BestbeforePoint
        self.chazhi_PointSelect[transSSTime:transSSTime + MinNumPoints, :] = chazhi_PointSelect.T
        self.transSSTime[TransTime]=transSSTime + MinNumPoints#每一段结束后总共共多少个点(统计所有圈数)
        self.tSSparaTime[self.TranslapTime, self.overlap]=self.tSSparaTime[self.TranslapTime-1, self.overlap]+ MinNumPoints  #每一段结束后，当前圈数总共有多少个点
        self.TranslapTime=self.TranslapTime+1 #TranslapTime表示这一圈的第几小段轨迹

        self.zVector=self.reallyroadPoint[10,:]
        #快到达终点处时处理
        if (self.zVector[4]  > map.TrackLength ): #为了制造圈与圈之间的起始点差，到终点时让经验轨迹倒退几步
           # houtui=int(np.floor(self.reallyroadPoint.shape[0]/3))  #到终点时后退1/3
           # self.zVector = self.reallyroadPoint[-houtui, :]- map.TrackLength
            self.zVector[4]=self.zVector[4]-map.TrackLength
            self.tSSlapsTime[self.overlap] = self.transSSTime[TransTime]  # 第overlap圈的时候，总共有多少个点（所有圈数）
            self.TranslapTime = 0
         #   self.passway=1    #防止后退完回到终点再次被后退
            print "迁移第%d圈" % (self.overlap)
            self.overlap = self.overlap + 1
       # if (self.zVector[4]  > map.TrackLength ) and (self.passway==1): #后退完再回到终点时
       #     self.zVector[4] = self.zVector[4] - map.TrackLength #满一圈时要减去一圈的路程值
       #     self.passway=0























    def addTrajectory(self, ClosedLoopData,LMPController,Laps):
        """update iteration index and construct SS, uSS and Qfun
        Arguments:
            ClosedLoopData: ClosedLoopData object
        """
        it = self.it
      #  print ClosedLoopData.SimTimeTot
        for it in range(2,Laps):
            self.TimeSS[it] = LMPController.TimeSS[it]  #过终点线仍然+1
            self.LapCounter[it] = ClosedLoopData.SimTimeTot[it] #一圈的仿真时间（不+1）
            self.SS[0:(self.TimeSS[it] + 1), :, it] = LMPController.SS[0:(self.TimeSS[it] + 1), :,it]
            self.SS_glob[0:(self.TimeSS[it] + 1), :, it] = LMPController.SS_glob[0:(self.TimeSS[it] + 1), :,it]
            self.uSS[0:self.TimeSS[it], :, it] = LMPController.uSS[0:(self.TimeSS[it]), :,it]
          #  print self.SS_glob[0:(self.TimeSS[it] + 1), 4, it]
           # print self.SS_glob[0:3, 4, it]
            for i in np.arange(1, self.SS_glob.shape[0]-2):
                if self.SS_glob[i + 2, 4, it]==10000:
                    break
                x = self.SS_glob[i-1:i + 2, 4, it]
                y = self.SS_glob[i-1:i + 2, 5, it]
              #  print x
             #   print y
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
    return -kappa
#求两点间弧长
def Arclength(bcur1,bcur0,x1,y1,x0,y0):
    bcur=abs(0.5*(bcur1+bcur0))
    if bcur<=1e-6:
        arclen=np.sqrt((x1-x0)**2+(y1-y0)**2)
    else:
        r=1/bcur
        d=np.sqrt((x1-x0)**2+(y1-y0)**2)
        rad=np.arcsin(0.5*d*bcur)
        arclen=2*rad*r
    return arclen