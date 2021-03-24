# coding=UTF-8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches

import pdb

def plotTrajectory(map, x, x_glob, u):
    Points = int(np.floor(10 * (map.PointAndTangent[-1, 3] + map.PointAndTangent[-1, 4])))
    Points1 = np.zeros((Points, 2))
    Points2 = np.zeros((Points, 2))
    Points0 = np.zeros((Points, 2))
    for i in range(0, int(Points)):
        Points1[i, :] = map.getGlobalPosition(i * 0.1, map.halfWidth)
        Points2[i, :] = map.getGlobalPosition(i * 0.1, -map.halfWidth)
        Points0[i, :] = map.getGlobalPosition(i * 0.1, 0)
    plt.figure('SS0_xy')
    plt.plot(map.PointAndTangent[:, 0], map.PointAndTangent[:, 1], 'o')
    plt.plot(Points0[:, 0], Points0[:, 1], '--')
    plt.plot(Points1[:, 0], Points1[:, 1], '-b')
    plt.plot(Points2[:, 0], Points2[:, 1], '-b')
    plt.plot(x_glob[:, 4], x_glob[:, 5], '-r')
    #plt.savefig("C:/Users/EE506/Desktop/fig/lab/SS0_xy.svg") #svg格式能在Visio中处理并变成emf（visio）专用格式
    plt.savefig("D:/study/python project/lmpc/src/lab/SS0_xy.png")
    plt.show()
    plt.figure('SS0_state')
    plt.subplot(711)
    plt.plot(x[:, 4], x[:, 0], '-o')
    plt.ylabel('vx')
    plt.subplot(712)
    plt.plot(x[:, 4], x[:, 1], '-o')
    plt.ylabel('vy')
    plt.subplot(713)
    plt.plot(x[:, 4], x[:, 2], '-o')
    plt.ylabel('wz')
    plt.subplot(714)
    plt.plot(x[:, 4], x[:, 3], '-o')
    plt.ylabel('epsi')
    plt.subplot(715)
    plt.plot(x[:, 4], x[:, 5], '-o')
    plt.ylabel('ey')
    plt.subplot(716)
    plt.plot(x[0:-1, 4], u[:, 0], '-o')
    plt.ylabel('steering')
    plt.subplot(717)
    plt.plot(x[0:-1, 4], u[:, 1], '-o')
    plt.ylabel('acc')
    plt.show()




#生成所有圈数都叠在一块的LMPC图片
def plotClosedLoopLMPC(LMPController, map):
    SS_glob = LMPController.SS_glob
    LapCounter  = LMPController.LapCounter
    SS      = LMPController.SS
    uSS     = LMPController.uSS

    TotNumberIt = LMPController.it
    Points = int(np.floor(10 * (map.PointAndTangent[-1, 3] + map.PointAndTangent[-1, 4])))
    Points1 = np.zeros((Points, 2))
    Points2 = np.zeros((Points, 2))
    Points0 = np.zeros((Points, 2))
    for i in range(0, int(Points)):
        Points1[i, :] = map.getGlobalPosition(i * 0.1, map.halfWidth)
        Points2[i, :] = map.getGlobalPosition(i * 0.1, -map.halfWidth)
        Points0[i, :] = map.getGlobalPosition(i * 0.1, 0)

    plt.figure(1)
    plt.plot(map.PointAndTangent[:, 0], map.PointAndTangent[:, 1], 'o')
    plt.plot(Points0[:, 0], Points0[:, 1], '--')
    plt.plot(Points1[:, 0], Points1[:, 1], '-b')
    plt.plot(Points2[:, 0], Points2[:, 1], '-b')

    for i in range(2, TotNumberIt):
        plt.plot(SS_glob[0:LapCounter[i], 4, i], SS_glob[0:LapCounter[i], 5, i], '-r')
    plt.savefig("D:/study/python project/lmpc/src/lab/lmpc_xy.png")
    plt.figure(2)
    plt.subplot(711)
    for i in range(2, TotNumberIt):
        plt.plot(SS[0:LapCounter[i], 4, i], SS[0:LapCounter[i], 0, i], '-o')
    plt.ylabel('vx')
    plt.subplot(712)
    for i in range(2, TotNumberIt):
        plt.plot(SS[0:LapCounter[i], 4, i], SS[0:LapCounter[i], 1, i], '-o')
    plt.ylabel('vy')
    plt.subplot(713)
    for i in range(2, TotNumberIt):
        plt.plot(SS[0:LapCounter[i], 4, i], SS[0:LapCounter[i], 2, i], '-o')
    plt.ylabel('wz')
    plt.subplot(714)
    for i in range(2, TotNumberIt):
        plt.plot(SS[0:LapCounter[i], 4, i], SS[0:LapCounter[i], 3, i], '-o')
    plt.ylabel('epsi')
    plt.subplot(715)
    for i in range(2, TotNumberIt):
        plt.plot(SS[0:LapCounter[i], 4, i], SS[0:LapCounter[i], 5, i], '-o')
    plt.ylabel('ey')
    plt.subplot(716)
    for i in range(2, TotNumberIt):
        plt.plot(uSS[0:LapCounter[i] - 1, 0, i], '-o')
    plt.ylabel('Steering')
    plt.subplot(717)
    for i in range(2, TotNumberIt):
        plt.plot(uSS[0:LapCounter[i] - 1, 1, i], '-o')
    plt.ylabel('Acc')




def animation_xy(map, LMPCOpenLoopData, LMPController, it):
    SS_glob = LMPController.SS_glob
    LapCounter = LMPController.LapCounter
    SS = LMPController.SS
    uSS = LMPController.uSS

    Points = int(np.floor(10 * (map.PointAndTangent[-1, 3] + map.PointAndTangent[-1, 4])))
    Points1 = np.zeros((Points, 2))
    Points2 = np.zeros((Points, 2))
    Points0 = np.zeros((Points, 2))

    for i in range(0, int(Points)):
        Points1[i, :] = map.getGlobalPosition(i * 0.1, map.halfWidth)
        Points2[i, :] = map.getGlobalPosition(i * 0.1, -map.halfWidth)
        Points0[i, :] = map.getGlobalPosition(i * 0.1, 0)

    plt.figure(200)
    plt.plot(map.PointAndTangent[:, 0], map.PointAndTangent[:, 1], 'o')
    plt.plot(Points0[:, 0], Points0[:, 1], '--')
    plt.plot(Points1[:, 0], Points1[:, 1], '-b')
    plt.plot(Points2[:, 0], Points2[:, 1], '-b')
    plt.plot(SS_glob[0:LapCounter[it], 4, it], SS_glob[0:LapCounter[it], 5, it], '-ok', label="Closed-loop trajectory",zorder=-1)

    ax = plt.axes()
    SSpoints_x = []; SSpoints_y = []
    xPred = []; yPred = []
    SSpoints, = ax.plot(SSpoints_x, SSpoints_y, 'sb', label="SS",zorder=0)
    line, = ax.plot(xPred, yPred, '-or', label="Predicted Trajectory",zorder=1)

    v = np.array([[ 1.,  1.],
                  [ 1., -1.],
                  [-1., -1.],
                  [-1.,  1.]])
    rec = patches.Polygon(v, alpha=0.7,closed=True, fc='r', ec='k',zorder=10)
    ax.add_patch(rec)

    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3)

    N = LMPController.N
    numSS_Points = LMPController.numSS_Points
    for i in range(0, int(LMPController.LapCounter[it])):

        xPred = np.zeros((N+1, 1)); yPred = np.zeros((N+1, 1))
        SSpoints_x = np.zeros((numSS_Points, 1)); SSpoints_y = np.zeros((numSS_Points, 1))

        for j in range(0, N+1):
            xPred[j,0], yPred[j,0]  = map.getGlobalPosition( LMPCOpenLoopData.PredictedStates[j, 4, i, it],
                                                             LMPCOpenLoopData.PredictedStates[j, 5, i, it] )

            if j == 0:
                x = SS_glob[i, 4, it]
                y = SS_glob[i, 5, it]
                psi = SS_glob[i, 3, it]
                l = 0.4; w = 0.2
                car_x = [ x + l * np.cos(psi) - w * np.sin(psi), x + l*np.cos(psi) + w * np.sin(psi),
                          x - l * np.cos(psi) + w * np.sin(psi), x - l * np.cos(psi) - w * np.sin(psi)]
                car_y = [ y + l * np.sin(psi) + w * np.cos(psi), y + l * np.sin(psi) - w * np.cos(psi),
                          y - l * np.sin(psi) - w * np.cos(psi), y - l * np.sin(psi) + w * np.cos(psi)]




        for j in range(0, numSS_Points):
            SSpoints_x[j,0], SSpoints_y[j,0] = map.getGlobalPosition(LMPCOpenLoopData.SSused[4, j, i, it],
                                                                     LMPCOpenLoopData.SSused[5, j, i, it])
        SSpoints.set_data(SSpoints_x, SSpoints_y)

        line.set_data(xPred, yPred)

        rec.set_xy(np.array([car_x, car_y]).T)

        plt.draw()
        plt.pause(1e-17)

def animation_states(map, LMPCOpenLoopData, LMPController, it):
    SS_glob = LMPController.SS_glob
    LapCounter = LMPController.LapCounter
    SS = LMPController.SS
    uSS = LMPController.uSS

    xdata = []; ydata = []
    fig = plt.figure(100)

    axvx = fig.add_subplot(3, 2, 1)
    plt.plot(SS[0:LapCounter[it], 4, it], SS[0:LapCounter[it], 0, it], '-ok', label="Closed-loop trajectory")
    lineSSvx, = axvx.plot(xdata, ydata, 'sb-', label="SS")
    linevx, = axvx.plot(xdata, ydata, 'or-', label="Predicted Trajectory")
    plt.ylabel("vx")
    plt.xlabel("s")

    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3)

    axvy = fig.add_subplot(3, 2, 2)
    axvy.plot(SS[0:LapCounter[it], 4, it], SS[0:LapCounter[it], 1, it], '-ok')
    lineSSvy, = axvy.plot(xdata, ydata, 'sb-')
    linevy, = axvy.plot(xdata, ydata, 'or-')
    plt.ylabel("vy")
    plt.xlabel("s")

    axwz = fig.add_subplot(3, 2, 3)
    axwz.plot(SS[0:LapCounter[it], 4, it], SS[0:LapCounter[it], 2, it], '-ok')
    lineSSwz, = axwz.plot(xdata, ydata, 'sb-')
    linewz, = axwz.plot(xdata, ydata, 'or-')
    plt.ylabel("wz")
    plt.xlabel("s")

    axepsi = fig.add_subplot(3, 2, 4)
    axepsi.plot(SS[0:LapCounter[it], 4, it], SS[0:LapCounter[it], 3, it], '-ok')
    lineSSepsi, = axepsi.plot(xdata, ydata, 'sb-')
    lineepsi, = axepsi.plot(xdata, ydata, 'or-')
    plt.ylabel("epsi")
    plt.xlabel("s")

    axey = fig.add_subplot(3, 2, 5)
    axey.plot(SS[0:LapCounter[it], 4, it], SS[0:LapCounter[it], 5, it], '-ok')
    lineSSey, = axey.plot(xdata, ydata, 'sb-')
    lineey, = axey.plot(xdata, ydata, 'or-')
    plt.ylabel("ey")
    plt.xlabel("s")

    Points = np.floor(10 * (map.PointAndTangent[-1, 3] + map.PointAndTangent[-1, 4])).astype(np.int)
    Points1 = np.zeros((Points, 2))
    Points2 = np.zeros((Points, 2))
    Points0 = np.zeros((Points, 2))
    for i in range(0, int(Points)):
        Points1[i, :] = map.getGlobalPosition(i * 0.1, map.halfWidth)
        Points2[i, :] = map.getGlobalPosition(i * 0.1, -map.halfWidth)
        Points0[i, :] = map.getGlobalPosition(i * 0.1, 0)

    axtr = fig.add_subplot(3, 2, 6)
    plt.plot(map.PointAndTangent[:, 0], map.PointAndTangent[:, 1], 'o')
    plt.plot(Points0[:, 0], Points0[:, 1], '--')
    plt.plot(Points1[:, 0], Points1[:, 1], '-b')
    plt.plot(Points2[:, 0], Points2[:, 1], '-b')
    plt.plot(SS_glob[0:LapCounter[it], 4, it], SS_glob[0:LapCounter[it], 5, it], '-ok')

    SSpoints_x = []; SSpoints_y = []
    xPred = []; yPred = []
    SSpoints_tr, = axtr.plot(SSpoints_x, SSpoints_y, 'sb')
    line_tr, = axtr.plot(xPred, yPred, '-or')

    N = LMPController.N
    numSS_Points = LMPController.numSS_Points
    for i in range(0, int(LMPController.LapCounter[it])):

        xPred    = LMPCOpenLoopData.PredictedStates[:, :, i, it]
        SSpoints = LMPCOpenLoopData.SSused[:, :, i, it]

        linevx.set_data(xPred[:, 4], xPred[:, 0]);   axvx.set_title(str(xPred[0, 0]))
        linevy.set_data(xPred[:, 4], xPred[:, 1]);   axvy.set_title(str(xPred[0, 1]))
        linewz.set_data(xPred[:, 4], xPred[:, 2]);   axwz.set_title(str(xPred[0, 2]))
        lineepsi.set_data(xPred[:, 4], xPred[:, 3]); axepsi.set_title(str(xPred[0, 3]))
        lineey.set_data(xPred[:, 4], xPred[:, 5]);   axey.set_title(str(xPred[0, 5]))

        epsiReal = xPred[0, 3]

        lineSSvx.set_data(SSpoints[4,:], SSpoints[0,:])
        lineSSvy.set_data(SSpoints[4,:], SSpoints[1,:])
        lineSSwz.set_data(SSpoints[4,:], SSpoints[2,:])
        lineSSepsi.set_data(SSpoints[4,:], SSpoints[3,:])
        lineSSey.set_data(SSpoints[4,:], SSpoints[5,:])

        xPred = np.zeros((N + 1, 1));yPred = np.zeros((N + 1, 1))
        SSpoints_x = np.zeros((numSS_Points, 1));SSpoints_y = np.zeros((numSS_Points, 1))

        for j in range(0, N + 1):
            xPred[j, 0], yPred[j, 0] = map.getGlobalPosition(LMPCOpenLoopData.PredictedStates[j, 4, i, it],
                                                             LMPCOpenLoopData.PredictedStates[j, 5, i, it])

        for j in range(0, numSS_Points):
            SSpoints_x[j, 0], SSpoints_y[j, 0] = map.getGlobalPosition(LMPCOpenLoopData.SSused[4, j, i, it],
                                                                       LMPCOpenLoopData.SSused[5, j, i, it])

        line_tr.set_data(xPred, yPred)


        vec = np.array([xPred[0, 0], yPred[0, 0]]) - np.array([SS_glob[i, 4, it], SS_glob[i, 5, it]])

        s, ey, epsi, _ = map.getLocalPosition( SS_glob[i, 4, it], SS_glob[i, 5, it], SS_glob[i, 3, it])
        axtr.set_title(str(s)+" "+str(ey)+" "+str(epsi))

        # axepsi.set_title(str(epsiReal)+" "+str(epsi))
        SSpoints_tr.set_data(SSpoints_x, SSpoints_y)

        plt.draw()
        plt.pause(1e-17)

def saveGif_xyResults(map, LMPCOpenLoopData, LMPController, it):
    SS_glob = LMPController.SS_glob
    LapCounter = LMPController.LapCounter
    SS = LMPController.SS
    uSS = LMPController.uSS

    Points = int(np.floor(10 * (map.PointAndTangent[-1, 3] + map.PointAndTangent[-1, 4])))  #道路的最后一个点扩大十倍并向下取整：如将3.21变成32
    Points1 = np.zeros((Points, 2))
    Points2 = np.zeros((Points, 2))
    Points0 = np.zeros((Points, 2))
    for i in range(0, int(Points)):
        Points1[i, :] = map.getGlobalPosition(i * 0.1, map.halfWidth)  #绘制道路的边界线
        Points2[i, :] = map.getGlobalPosition(i * 0.1, -map.halfWidth) #绘制道路的另一条边界线
        Points0[i, :] = map.getGlobalPosition(i * 0.1, 0)   #绘制道路的中心线

    fig = plt.figure(101)
    # plt.ylim((-5, 1.5))
    fig.set_tight_layout(True)  #自动紧凑布局
    plt.plot(map.PointAndTangent[:, 0], map.PointAndTangent[:, 1], 'o')
    #蓝色的路宽两条线及路中间的虚线
    plt.plot(Points0[:, 0], Points0[:, 1], '--')
    plt.plot(Points1[:, 0], Points1[:, 1], '-b')
    plt.plot(Points2[:, 0], Points2[:, 1], '-b')
    #黑色的闭环轨迹
    plt.plot(SS_glob[0:LapCounter[it], 4, it], SS_glob[0:LapCounter[it], 5, it], '-ok', label="Closed-loop trajectory", markersize=1,zorder=-1)

    ax = plt.axes()
    SSpoints_x = []; SSpoints_y = []
    xPred = []; yPred = []
    SSpoints, = ax.plot(SSpoints_x, SSpoints_y, 'og', label="SS",zorder=0)   #ss集绘制
    line, = ax.plot(xPred, yPred, '-or', label="Predicted Trajectory",zorder=1)  #预测的状态变量绘制

    v = np.array([[ 1.,  1.],
                  [ 1., -1.],
                  [-1., -1.],
                  [-1.,  1.]])
    #绿色的小车模型
    rec = patches.Polygon(v, alpha=0.7,closed=True, fc='w', ec='k',zorder=10)
    ax.add_patch(rec)

    plt.legend(mode="expand", ncol=3)
    # plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
    #             mode="expand", borderaxespad=0, ncol=3)

    N = LMPController.N
    numSS_Points = LMPController.numSS_Points

    def update(i):
        xPred = np.zeros((N + 1, 1)); yPred = np.zeros((N + 1, 1))
        SSpoints_x = np.zeros((numSS_Points, 1)); SSpoints_y = np.zeros((numSS_Points, 1))
        #预测的状态变量
        for j in range(0, N + 1):
            xPred[j, 0], yPred[j, 0] = map.getGlobalPosition(LMPCOpenLoopData.PredictedStates[j, 4, i, it],
                                                             LMPCOpenLoopData.PredictedStates[j, 5, i, it])

            if j == 0:
                x = SS_glob[i, 4, it]
                y = SS_glob[i, 5, it]
                psi = SS_glob[i, 3, it]
                l = 0.4;w = 0.2
                car_x = [x + l * np.cos(psi) - w * np.sin(psi), x + l * np.cos(psi) + w * np.sin(psi),
                         x - l * np.cos(psi) + w * np.sin(psi), x - l * np.cos(psi) - w * np.sin(psi)]
                car_y = [y + l * np.sin(psi) + w * np.cos(psi), y + l * np.sin(psi) - w * np.cos(psi),
                         y - l * np.sin(psi) - w * np.cos(psi), y - l * np.sin(psi) + w * np.cos(psi)]
        #ss集
        for j in range(0, numSS_Points):
            SSpoints_x[j, 0], SSpoints_y[j, 0] = map.getGlobalPosition(LMPCOpenLoopData.SSused[4, j, i, it],
                                                                       LMPCOpenLoopData.SSused[5, j, i, it])
        SSpoints.set_data(SSpoints_x, SSpoints_y)

        line.set_data(xPred, yPred)

        rec.set_xy(np.array([car_x, car_y]).T)

    anim = FuncAnimation(fig, update, frames=np.arange(0, int(LMPController.LapCounter[it])), interval=100)

    anim.save('ClosedLoop.mp4', dpi=80, writer='imagemagick')



def Save_statesAnimation(map, LMPCOpenLoopData, LMPController, it):
    SS_glob = LMPController.SS_glob
    LapCounter = LMPController.LapCounter
    SS = LMPController.SS
    uSS = LMPController.uSS

    xdata = []; ydata = []
    fig = plt.figure()
    fig.set_tight_layout(True)

    axvx = fig.add_subplot(3, 2, 1)
    plt.plot(SS[0:LapCounter[it], 4, it], SS[0:LapCounter[it], 0, it], '-ok', label="Closed-loop trajectory")
    lineSSvx, = axvx.plot(xdata, ydata, 'sb-', label="SS")
    linevx, = axvx.plot(xdata, ydata, 'or-', label="Predicted Trajectory")
    plt.ylabel("vx")
    plt.xlabel("s")

    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3)

    axvy = fig.add_subplot(3, 2, 2)
    axvy.plot(SS[0:LapCounter[it], 4, it], SS[0:LapCounter[it], 1, it], '-ok')
    lineSSvy, = axvy.plot(xdata, ydata, 'sb-')
    linevy, = axvy.plot(xdata, ydata, 'or-')
    plt.ylabel("vy")
    plt.xlabel("s")

    axwz = fig.add_subplot(3, 2, 3)
    axwz.plot(SS[0:LapCounter[it], 4, it], SS[0:LapCounter[it], 2, it], '-ok')
    lineSSwz, = axwz.plot(xdata, ydata, 'sb-')
    linewz, = axwz.plot(xdata, ydata, 'or-')
    plt.ylabel("wz")
    plt.xlabel("s")

    axepsi = fig.add_subplot(3, 2, 4)
    axepsi.plot(SS[0:LapCounter[it], 4, it], SS[0:LapCounter[it], 3, it], '-ok')
    lineSSepsi, = axepsi.plot(xdata, ydata, 'sb-')
    lineepsi, = axepsi.plot(xdata, ydata, 'or-')
    plt.ylabel("epsi")
    plt.xlabel("s")

    axey = fig.add_subplot(3, 2, 5)
    axey.plot(SS[0:LapCounter[it], 4, it], SS[0:LapCounter[it], 5, it], '-ok')
    lineSSey, = axey.plot(xdata, ydata, 'sb-')
    lineey, = axey.plot(xdata, ydata, 'or-')
    plt.ylabel("ey")
    plt.xlabel("s")

    Points = int(np.floor(10 * (map.PointAndTangent[-1, 3] + map.PointAndTangent[-1, 4])))
    Points1 = np.zeros((Points, 2))
    Points2 = np.zeros((Points, 2))
    Points0 = np.zeros((Points, 2))
    for i in range(0, int(Points)):
        Points1[i, :] = map.getGlobalPosition(i * 0.1, map.halfWidth)
        Points2[i, :] = map.getGlobalPosition(i * 0.1, -map.halfWidth)
        Points0[i, :] = map.getGlobalPosition(i * 0.1, 0)

    axtr = fig.add_subplot(3, 2, 6)
    plt.plot(map.PointAndTangent[:, 0], map.PointAndTangent[:, 1], 'o')
    plt.plot(Points0[:, 0], Points0[:, 1], '--')
    plt.plot(Points1[:, 0], Points1[:, 1], '-b')
    plt.plot(Points2[:, 0], Points2[:, 1], '-b')
    plt.plot(SS_glob[0:LapCounter[it], 4, it], SS_glob[0:LapCounter[it], 5, it], '-ok')

    SSpoints_x = []; SSpoints_y = []
    xPred = []; yPred = []
    SSpoints_tr, = axtr.plot(SSpoints_x, SSpoints_y, 'sb')
    line_tr, = axtr.plot(xPred, yPred, '-or')

    N = LMPController.N
    numSS_Points = LMPController.numSS_Points

    def update(i):
        xPred    = LMPCOpenLoopData.PredictedStates[:, :, i, it]
        SSpoints = LMPCOpenLoopData.SSused[:, :, i, it]

        linevx.set_data(xPred[:, 4], xPred[:, 0])
        linevy.set_data(xPred[:, 4], xPred[:, 1])
        linewz.set_data(xPred[:, 4], xPred[:, 2])
        lineepsi.set_data(xPred[:, 4], xPred[:, 3])
        lineey.set_data(xPred[:, 4], xPred[:, 5])

        lineSSvx.set_data(SSpoints[4,:], SSpoints[0,:])
        lineSSvy.set_data(SSpoints[4,:], SSpoints[1,:])
        lineSSwz.set_data(SSpoints[4,:], SSpoints[2,:])
        lineSSepsi.set_data(SSpoints[4,:], SSpoints[3,:])
        lineSSey.set_data(SSpoints[4,:], SSpoints[5,:])

        xPred = np.zeros((N + 1, 1));yPred = np.zeros((N + 1, 1))
        SSpoints_x = np.zeros((numSS_Points, 1));SSpoints_y = np.zeros((numSS_Points, 1))

        for j in range(0, N + 1):
            xPred[j, 0], yPred[j, 0] = map.getGlobalPosition(LMPCOpenLoopData.PredictedStates[j, 4, i, it],
                                                             LMPCOpenLoopData.PredictedStates[j, 5, i, it])

        for j in range(0, numSS_Points):
            SSpoints_x[j, 0], SSpoints_y[j, 0] = map.getGlobalPosition(LMPCOpenLoopData.SSused[4, j, i, it],
                                                                       LMPCOpenLoopData.SSused[5, j, i, it])

        line_tr.set_data(xPred, yPred)
        SSpoints_tr.set_data(SSpoints_x, SSpoints_y)

    anim = FuncAnimation(fig, update, frames=np.arange(0, int(LMPController.LapCounter[it])), interval=100)

    anim.save('ClosedLoopStates.mp4', dpi=80, writer='imagemagick')

def TLMPC_xyResults(map, TLMPCOpenLoopData, TLMPController):
   # SS_glob = LMPController.SS_glob
    LapCounter = TLMPController.LapCounter
    TransTime= TLMPController.TransTime+1
  #  SS = LMPController.SS
  #  uSS = LMPController.uSS
#绘制道路
    Points = int(np.floor(10 * (map.PointAndTangent[-1, 3] + map.PointAndTangent[-1, 4])))  #道路的最后一个点扩大十倍并向下取整：如将3.21变成32
    Points1 = np.zeros((Points, 2))
    Points2 = np.zeros((Points, 2))
    Points0 = np.zeros((Points, 2))
    for i in range(0, int(Points)):
        Points1[i, :] = map.getGlobalPosition(i * 0.1, map.halfWidth)  #绘制道路的边界线
        Points2[i, :] = map.getGlobalPosition(i * 0.1, -map.halfWidth) #绘制道路的另一条边界线
        Points0[i, :] = map.getGlobalPosition(i * 0.1, 0)   #绘制道路的中心线

    fig = plt.figure(101)
    # plt.ylim((-5, 1.5))
    fig.set_tight_layout(True)  #自动紧凑布局
    plt.plot(map.PointAndTangent[:, 0], map.PointAndTangent[:, 1], 'o')
    #蓝色的路宽两条线及路中间的虚线
    plt.plot(Points0[:, 0], Points0[:, 1], '--')
    plt.plot(Points1[:, 0], Points1[:, 1], '-b')
    plt.plot(Points2[:, 0], Points2[:, 1], '-b')
#绘制运动轨迹bel="Closed-loop trajectory", markersize=1,zorder=-1)

   #绘制静态sr集

    numSR_Points = TLMPController.numSR_points  #总点数
    #黑色的闭环轨迹
   # plt.plot(SS_glob[0:LapCounter[it], 4, it], SS_glob[0:LapCounter[it], 5, it], '-ok', la
    SRpoints_x = np.empty((numSR_Points,TransTime));SRpoints_y =  np.empty((numSR_Points,TransTime)) #相似路径
    RRpoints_x = np.empty((numSR_Points,TransTime));RRpoints_y = np.empty((numSR_Points,TransTime))   #实际道路中心线
  #  SRpoints1_x = np.empty(numSR_Points);SRpoints1_y = np.empty(numSR_Points)  #相似路径虚拟边界1
   # SRpoints2_x = np.empty(numSR_Point                                          s);SRpoints2_y =  np.empty(numSR_Points)  #相似路径虚拟边界2
    for i in range(0,TransTime):
      for j in range(0, numSR_Points):
    #     print 'ceshi'
         if TLMPCOpenLoopData.RRused[4, j, i]<0:
             TLMPCOpenLoopData.RRused[4,j,i]+=map.TrackLength
        #选出的相似路径
         SRpoints_x[j,i], SRpoints_y[j,i] = map.getGlobalPosition(TLMPCOpenLoopData.TSRused[4, j, i],
                                                               TLMPCOpenLoopData.TSRused[5, j, i])
   #     #相似路径的虚拟路径边界
   #      SRpoints1_x[j], SRpoints1_y[j] = map.getGlobalPosition(TLMPCOpenLoopData.TSRused[4, j, 0, it],
         #                                            TLMPCOpenLoopData.TSRused[5, j, 0, it]+map.halfWidth)
     #    SRpoints2_x[j], SRpoints2_y[j] = map.getGlobalPosition(TLMPCOpenLoopData.TSRused[4, j, 0, it],
        #                                                       TLMPCOpenLoopData.TSRused[5, j, 0, it]-map.halfWidth)
        #实际路径中心线
         RRpoints_x[j,i], RRpoints_y[j,i] = map.getGlobalPosition(TLMPCOpenLoopData.RRused[4, j, i],
                                                              TLMPCOpenLoopData.RRused[5, j, i])
      if i==0:
      #   plt.plot(SRpoints_x[:,i], SRpoints_y[:,i],'--r',label="Virtual centerline")
         plt.plot(RRpoints_x[-1,i], RRpoints_y[-1,i], 'o',color='red',label="current centerline")
         plt.plot(RRpoints_x[0, i], RRpoints_y[0, i], 'o', color='red')
      else:
      #   plt.plot(SRpoints_x[:, i], SRpoints_y[:, i], '--r')
         plt.plot(RRpoints_x[-1, i], RRpoints_y[-1, i], 'o',color='red')
    #  if i==2:
    #      plt.plot(SRpoints_x[:, i], SRpoints_y[:, i], '--r')

      if i==0:
          transSSTime=TLMPController.transSSTime[i]
      else:
          transSSTime = TLMPController.transSSTime[i]-TLMPController.transSSTime[i-1]

      tSSpoints_x = np.empty(transSSTime);
      tSSpoints_y = np.empty(transSSTime)  # 迁移后的最优轨迹
      BeforeSSpoints_x = np.empty(transSSTime);
      BeforeSSpoints_y = np.empty(transSSTime)  # 迁移后的最优轨迹
      chazhiSSpoints_x = np.empty(transSSTime);
      chazhiSSpoints_y = np.empty(transSSTime)
      for j in range(0, transSSTime):
         if i==0:
             jj=j
         else:
             jj=TLMPController.transSSTime[i-1]+j  #算上之前路径段后j的位置
         #迁移后的最优轨迹
         tSSpoints_x[j], tSSpoints_y[j] = map.getGlobalPosition(TLMPController.transSS[jj,4],
                                                     TLMPController.transSS[jj,5])
         BeforeSSpoints_x[j], BeforeSSpoints_y[j] = map.getGlobalPosition(TLMPController.BestbeforePoint[jj, 4],
                                                       TLMPController.BestbeforePoint[jj, 5])
         chazhiSSpoints_x[j], chazhiSSpoints_y[j] = map.getGlobalPosition(TLMPController.chazhi_PointSelect[jj, 4],
                                                       TLMPController.chazhi_PointSelect[jj, 5])
      if i==0:
        plt.plot(tSSpoints_x, tSSpoints_y, linestyle='-', color='green', label="Transfered Track")
      #  plt.plot(BeforeSSpoints_x, BeforeSSpoints_y, linestyle='--', color='green', label="Before Track")
      else:
        plt.plot(tSSpoints_x, tSSpoints_y, linestyle='-', color='green',)
    #  if i==2:
    #    plt.plot(BeforeSSpoints_x, BeforeSSpoints_y, linestyle='--', color='green',)
   # plt.plot(BeforeSSpoints_x, BeforeSSpoints_y, linestyle='--', color='green', label="Before Track") 此注释为例子

 #   plt.plot(SRpoints1_x, SRpoints1_y, '--b',label="Virtual boundaries")
 #   plt.plot(SRpoints2_x, SRpoints2_y, '--b')

  #  plt.plot(chazhiSSpoints_x, chazhiSSpoints_y, '-ob', label="chazhi")
    plt.legend(loc='best')

def newsaveGif_xyResults(map, LMPCOpenLoopData, LMPController, it,LMPCloseloopData):

    SS_glob = LMPController.SS_glob
    LapCounter = LMPController.LapCounter
    x_glob=LMPCloseloopData.x_glob
    SS = LMPController.SS
    uSS = LMPController.uSS

    Points = int(np.floor(10 * (map.PointAndTangent[-1, 3] + map.PointAndTangent[-1, 4])))  #道路的最后一个点扩大十倍并向下取整：如将3.21变成32
    Points1 = np.zeros((Points, 2))
    Points2 = np.zeros((Points, 2))
    Points0 = np.zeros((Points, 2))

    for i in range(0, int(Points)):
        Points1[i, :] = map.getGlobalPosition(i * 0.1, map.halfWidth)  #绘制道路的边界线
        Points2[i, :] = map.getGlobalPosition(i * 0.1, -map.halfWidth) #绘制道路的另一条边界线
        Points0[i, :] = map.getGlobalPosition(i * 0.1, 0)   #绘制道路的中心线

    fig = plt.figure(101)
    # plt.ylim((-5, 1.5))
    fig.set_tight_layout(True)  #自动紧凑布局
    plt.plot(map.PointAndTangent[:, 0], map.PointAndTangent[:, 1], 'o')
    #蓝色的路宽两条线及路中间的虚线
    plt.plot(Points0[:, 0], Points0[:, 1], '--')
    plt.plot(Points1[:, 0], Points1[:, 1], '-b')
    plt.plot(Points2[:, 0], Points2[:, 1], '-b')
    #黑色的闭环轨迹
    plt.plot(x_glob[:, 4], x_glob[:, 5], '-ok', label="Closed-loop trajectory", markersize=1,zorder=-1)
    ax = plt.axes()

    SSpoints_x = []; SSpoints_y = []
    xPred = []; yPred = []
    SSpoints, = ax.plot(SSpoints_x, SSpoints_y, 'og', label="SS",zorder=0)   #ss集绘制
    line, = ax.plot(xPred, yPred, '-or', label="Predicted Trajectory",zorder=1)  #预测的状态变量绘制

    v = np.array([[ 1.,  1.],
                  [ 1., -1.],
                  [-1., -1.],
                  [-1.,  1.]])

    plt.legend(mode="expand", ncol=3)
    # plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
    #             mode="expand", borderaxespad=0, ncol=3)

    N = LMPController.N
    numSS_Points = LMPController.numSS_Points


    def update(i):


        xPred = np.zeros((N + 1, 1)); yPred = np.zeros((N + 1, 1))
        SSpoints_x = np.zeros((numSS_Points, 1)); SSpoints_y = np.zeros((numSS_Points, 1))
        #预测的状态变量
        for j in range(0, N + 1):

             xPred[j, 0], yPred[j, 0] = map.getGlobalPosition(LMPCOpenLoopData.PredictedStates[j, 4, i, it],
                                                             LMPCOpenLoopData.PredictedStates[j, 5, i, it])


        #ss集
        for j in range(0, numSS_Points):

                SSpoints_x[j, 0], SSpoints_y[j, 0] = map.getGlobalPosition(LMPCOpenLoopData.SSused[4, j, i, it],
                                                                       LMPCOpenLoopData.SSused[5, j, i, it])
        SSpoints.set_data(SSpoints_x, SSpoints_y)

        line.set_data(xPred, yPred)


    anim = FuncAnimation(fig, update, frames=np.arange(0, 200), interval=100)

    anim.save('ClosedLoop.mp4', dpi=80, writer='imagemagick')