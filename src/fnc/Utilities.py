# coding=UTF-8
import numpy as np
import pdb
def Regression(x, u, lamb):
    """Estimates linear system dynamics
    x, u: date used in the regression
    lamb: regularization coefficient
    """

    # Want to solve W^* = argmin sum_i ||W^T z_i - y_i ||_2^2 + lamb ||W||_F,
    # with z_i = [x_i u_i] and W \in R^{n + d} x n
    Y = x[2:x.shape[0], :]                                             #输出的状态变量
    X = np.hstack((x[1:(x.shape[0] - 1), :], u[1:(x.shape[0] - 1), :])) #hatack将元素横着堆叠，输入的状态变量和输入u堆叠，用来计算输出

    Q = np.linalg.inv(np.dot(X.T, X) + lamb * np.eye(X.shape[1]))
    b = np.dot(X.T, Y)
    W = np.dot(Q, b)

    A = W.T[:, 0:6]   #状态变量
    B = W.T[:, 6:8]   #输出变量

    ErrorMatrix = np.dot(X, W) - Y   #X.T*Y*X/X.T*X=Y
    ErrorMax = np.max(ErrorMatrix, axis=0)
    ErrorMin = np.min(ErrorMatrix, axis=0)
    Error = np.vstack((ErrorMax, ErrorMin)) #最大误差和最小误差竖着堆叠

    return A, B, Error

def Curvature(s, PointAndTangent):
    """curvature computation
    s: curvilinear abscissa at which the curvature has to be evaluated
    PointAndTangent: points and tangent vectors defining the map (these quantities are initialized in the map object)
    """
    TrackLength = PointAndTangent[-1,3]+PointAndTangent[-1,4] #pointandtangent[ ,3]和[ ,4]分别表示cumulative s 和 segment length

    # In case on a lap after the first one
    while (s > TrackLength):
        s = s - TrackLength

    # Given s \in [0, TrackLength] compute the curvature
    # Compute the segment in which system is evolving
    #np.all进行交运算，找到s在哪个片段中
    index = np.all([[s >= PointAndTangent[:, 3]], [s < PointAndTangent[:, 3] + PointAndTangent[:, 4]]], axis=0)
   # print index
    #找到该片段的位置
    i = int(np.where(np.squeeze(index))[0])
   # print i
    curvature = PointAndTangent[i, 5]  #第一个s值的曲率

    return curvature

def getAngle(s, epsi, PointAndTangent):
    """curvature computation
    s: curvilinear abscissa at which the curvature has to be evaluated
    PointAndTangent: points and tangent vectors defining the map (these quantities are initialized in the map object)
    """
    TrackLength = PointAndTangent[-1,3]+PointAndTangent[-1,4]

    # In case on a lap after the first one
    while (s > TrackLength):
        s = s - TrackLength

    # Given s \in [0, TrackLength] compute the curvature
    # Compute the segment in which system is evolving
    index = np.all([[s >= PointAndTangent[:, 3]], [s < PointAndTangent[:, 3] + PointAndTangent[:, 4]]], axis=0)
    i = int(np.where(np.squeeze(index))[0])

    if i > 0:
        ang = PointAndTangent[i - 1, 2]
    else:
        ang = 0

    if PointAndTangent[i, 5] == 0:
        r= 0
    else:
        r = 1 / PointAndTangent[i, 5]  # Radius of curvature

    if r == 0:
        # On a straight part of the circuit
        angle_at_s = ang + epsi
    else:
        # On a curve
        cumulative_s = PointAndTangent[i, 3]
        relative_s = s - cumulative_s
        spanAng = relative_s / np.abs(r)  # Angle spanned by the circle
        psi = wrap(ang + spanAng * np.sign(r))  # Angle of the tangent vector at the last point of the segment
        # pdb.set_trace()
        angle_at_s = psi + epsi

    return angle_at_s

def wrap(angle):
    if angle < -np.pi:
        w_angle = 2 * np.pi + angle
    elif angle > np.pi:
        w_angle = angle - 2 * np.pi
    else:
        w_angle = angle

    return w_angle
