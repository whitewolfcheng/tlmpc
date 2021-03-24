#coding=utf-8
import time

class VirtualSimulator():
    def __init__(self, map, lap = 0, flagTLMPC=0):
        self.map = map
        self.laps = lap
        self.flagTLMPC=flagTLMPC

    def Sim(self, ClosedLoopData, Controller, TLMPCprediction=0,Tlaps=0):

        x      = ClosedLoopData.x
        x_glob = ClosedLoopData.x_glob
        u      = ClosedLoopData.u

        SimulationTime = 0

        for i in range(0, int(ClosedLoopData.Points)):  #在[0,points)中遍历， self.Points = int(Time / dt)，表示仿真中的点个数p0

            Controller.TransTime = i
            ceshistarttime=time.time()
            Controller.solve(x[i, :])                   #x中第一个量为预测个数，第二个量为6个状态变量
            ceshiendtime = time.time()

            TLMPCprediction.TSRused[:, :, i ]         = Controller.TSR_PointSelectedTot.T
            TLMPCprediction.RRused[:,:,i]             =Controller.reallyroadPoint.T

            if Controller.overflag==1:
                break
            if Controller.overlap>=Tlaps:
                break

        ClosedLoopData.SimTime = SimulationTime
        print "Number of laps completed: ", Controller.overlap#floor为向下取整，比如floor（2.1）=2.-1表示取x的最后一次预测