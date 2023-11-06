import math
import numpy as np
import matplotlib.pyplot as plt

# ISA数据计算
def ISA_calc(h):
    H_lowest, H0, H1, H2, H3, H4, H5, H6, H7 = -610, 0, 11000, 20000, 32000, 47000, 51000, 71000, 84852
    R, M, Re, Gamma = 8.3144598, 0.0289644, 6371000, 1.4
    R_special = R / M
    T0, P0, Rho0, g0 = 288.15, 1013.25, 1.225, 9.80665

    Lambda0to1, Lambda2to3, Lambda3to4, Lambda5to6, Lambda6to7 = -0.0065, 0.0010, 0.0028, -0.0028, -0.0020

    T1to2 = T0 + Lambda0to1 * H1
    T3 = T1to2 + Lambda2to3 * (H3 - H2)
    T4to5 = T3 + Lambda3to4 * (H4 - H3)
    T6 = T4to5 + Lambda5to6 * (H6 - H5)

    P1 = P0 * (1 + Lambda0to1 * (H1 - H0) / T0)**(-g0 / (R_special * Lambda0to1))
    P2 = P1 * math.exp(-g0 * (H2 - H1) / (R_special * T1to2))
    P3 = P2 * (1 + Lambda2to3 * (H3 - H2) / T1to2)**(-g0 / (R_special * Lambda2to3))
    P4 = P3 * (1 + Lambda3to4 * (H4 - H3) / T3)**(-g0 / (R_special * Lambda3to4))
    P5 = P4 * math.exp(-g0 * (H5 - H4) / (R_special * T4to5))
    P6 = P5 * (1 + Lambda5to6 * (H6 - H5) / T4to5)**(-g0 / (R_special * Lambda5to6))

    Rho1 = Rho0 * (1 + Lambda0to1 * (H1 - H0) / T0)**(-g0 / (R_special * Lambda0to1) - 1)
    Rho2 = Rho1 * math.exp(-g0 * (H2 - H1) / (R_special * T1to2))
    Rho3 = Rho2 * (1 + Lambda2to3 * (H3 - H2) / T1to2)**(-g0 / (R_special * Lambda2to3) - 1)
    Rho4 = Rho3 * (1 + Lambda3to4 * (H4 - H3) / T3)**(-g0 / (R_special * Lambda3to4) - 1)
    Rho5 = Rho4 * math.exp(-g0 * (H5 - H4) / (R_special * T4to5))
    Rho6 = Rho5 * (1 + Lambda5to6 * (H6 - H5) / T4to5)**(-g0 / (R_special * Lambda5to6) - 1)

    if h < H_lowest:
        print()
        print('Error: The altitude should between -610m ~ 51000m')
        print('Your altitude is', h, 'm')
        print('T P Rho a will be -1')
        print()
        T, P, Rho, a, g = -1, -1, -1, -1, -1
    elif h < H0:
        print()
        print('Warning0: The altitude is lower than sea level')
        print('Your altitude is', h, 'm')
        print()
        T = T0 + Lambda0to1 * (h - H0)
        P = P0 * (1 + Lambda0to1 * (h - H0) / T0)**(-g0 / (R_special * Lambda0to1))
        Rho = Rho0 * (1 + Lambda0to1 * (h - H0) / T0)**(-g0 / (R_special * Lambda0to1) - 1)
        a = math.sqrt(Gamma * R_special * T)
        g = g0 * (Re / (Re + h))**2
    elif h < H1:
        T = T0 + Lambda0to1 * (h - H0)
        P = P0 * (1 + Lambda0to1 * (h - H0) / T0)**(-g0 / (R_special * Lambda0to1))
        Rho = Rho0 * (1 + Lambda0to1 * (h - H0) / T0)**(-g0 / (R_special * Lambda0to1) - 1)
        a = math.sqrt(Gamma * R_special * T)
        g = g0 * (Re / (Re + h))**2
    elif h < H2:
        T = T1to2
        P = P1 * math.exp(-g0 * (h - H1) / (R_special * T1to2))
        Rho = Rho1 * math.exp(-g0 * (h - H1) / (R_special * T1to2))
        a = math.sqrt(Gamma * R_special * T)
        g = g0 * (Re / (Re + h))**2
    elif h < H3:
        T = T1to2 + Lambda2to3 * (h - H2)
        P = P2 * (1 + Lambda2to3 * (h - H2) / T1to2)**(-g0 / (R_special * Lambda2to3))
        Rho = Rho2 * (1 + Lambda2to3 * (h - H2) / T1to2)**(-g0 / (R_special * Lambda2to3) - 1)
        a = math.sqrt(Gamma * R_special * T)
        g = g0 * (Re / (Re + h))**2
    elif h < H4:
        T = T3 + Lambda3to4 * (h - H3)
        P = P3 * (1 + Lambda3to4 * (h - H3) / T3)**(-g0 / (R_special * Lambda3to4))
        Rho = Rho3 * (1 + Lambda3to4 * (h - H3) / T3)**(-g0 / (R_special * Lambda3to4) - 1)
        a = math.sqrt(Gamma * R_special * T)
        g = g0 * (Re / (Re + h))**2
    elif h < H5:
        T = T4to5
        P = P4 * math.exp(-g0 * (h - H4) / (R_special * T4to5))
        Rho = Rho4 * math.exp(-g0 * (h - H4) / (R_special * T4to5))
        a = math.sqrt(Gamma * R_special * T)
        g = g0 * (Re / (Re + h))**2
    elif h < H6:
        T = T4to5 + Lambda5to6 * (h - H5)
        P = P5 * (1 + Lambda5to6 * (h - H5) / T4to5)**(-g0 / (R_special * Lambda5to6))
        Rho = Rho5 * (1 + Lambda5to6 * (h - H5) / T4to5)**(-g0 / (R_special * Lambda5to6) - 1)
        a = math.sqrt(Gamma * R_special * T)
        g = g0 * (Re / (Re + h))**2
    elif h <= H7:
        T = T6 + Lambda6to7 * (h - H6)
        P = P6 * (1 + Lambda6to7 * (h - H6) / T6)**(-g0 / (R_special * Lambda6to7))
        Rho = Rho6 * (1 + Lambda6to7 * (h - H6) / T6)**(-g0 / (R_special * Lambda6to7) - 1)
        a = math.sqrt(Gamma * R_special * T)
        g = g0 * (Re / (Re + h))**2
    else:
        print()
        print('Error: The altitude should between -610m ~ 51000m')
        print('Your altitude is', h, 'm')
        print('T P Rho a will be -2')
        print()
        T, P, Rho, a, g = -2, -2, -2, -2, -2

    return T, P, Rho, a, g


'''
ALT: 给定高度 Altitude ( 英尺 | ft )
W: 飞机重量 Airplane Weight ( 磅 | lbs )
S: 机翼面积 Wing Area ( 平方英尺 | ft² )
CD0: 零升阻力系数 ( 无量纲 )
Epslion: 升致阻力因子
Eta: 螺旋桨效率
Psmax: 最大功率 ( 英尺·磅/秒 | ft·lbs/s )
m: 无量纲空气密度指数 Dimensionless Air Density Exponent ( 无量纲 )
CLmax: 最大升力系数 ( 无量纲 )

V: 给定的速度 ( 英尺/秒 | ft/s )
PWR: 给定的速度、高度和飞行条件下的平飞所需功率 ( 英尺·磅/秒 | ft·lbs/s )
PWRmax: 给定的高度和飞行条件下的最大功率 ( 英尺·磅/秒 | ft·lbs/s )

Rhos: ISA海平面 大气密度 ( 千克/立方米 | kg/m³ )
Rho: ISA模型下 给定高度 的 大气密度 ( 千克/立方米 | kg/m³ )

Vs: 失速速度 ( 英尺/秒 | ft/s )
Vmn: 最小需用功率速度(平飞有利速度) ( 英尺/秒 | ft/s )
Vmin: 最小平飞速度 ( 英尺/秒 | ft/s )
Vmax: 最大平飞速度 ( 英尺/秒 | ft/s )
Ts Te Tmin Tmax: 对应速度下平飞所需功率 ( 英尺·磅/秒 | ft·lbs/s )
'''

#数据准备
ALT, W, S, CD0, Epslion, Eta = 10000, 2900, 175, 0.026, 0.054, 0.8
Psmax, m, CLmax = 290*550, 0.6, 2.4
V = np.arange(20, 360, 10, float)
PWR = np.zeros(V.size)
PWRmax = np.zeros(V.size)

#数据计算
_, _, Rhos, _, _ = ISA_calc(0)
_, _, Rho, _, _ = ISA_calc(ALT * 0.3048)
Rhos*=0.0019
Rho*=0.0019

Vs = math.sqrt(2 * W / (Rho * S * CLmax))
Vmn = math.sqrt(2 * W / (Rho * S) * math.sqrt(Epslion / 3 / CD0))
roots = np.roots([1/2 * Rho * S * CD0, 0, 0, -Eta * Psmax * (Rho / Rhos)**m, 2 * Epslion * W**2 / Rho / S])
r = sorted(x.real for x in roots)
Vn1, Vn2, Vmin, Vmax = r
Tvs = 1/2 * Rho * Vs**3 * S * CD0 + 2 * Epslion * W**2 / (Rho * Vs * S)
Tvmn = 1/2 * Rho * Vmn**3 * S * CD0 + 2 * Epslion * W**2 / (Rho * Vmn * S)
Tvmin = 1/2 * Rho * Vmin**3 * S * CD0 + 2 * Epslion * W**2 / (Rho * Vmin * S)
Tvmax = 1/2 * Rho * Vmax**3 * S * CD0 + 2 * Epslion * W**2 / (Rho * Vmax * S)
    
print('Vs=%f, Vmn=%f, Vmin=%f, Vmax=%f' %(Vs, Vmn, Vmin, Vmax))

i = 0
for V_i in V:
    PWR[i]=1/2 * Rho * V_i**3 * S * CD0 + 2 * Epslion * W**2 / (Rho * V_i * S)
    PWRmax[i] = Eta * ((Rho / Rhos)**m) * Psmax
    i+=1

#画图
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(V, PWR)
ax.plot(V, PWRmax)
ax.scatter(Vs, Tvs, s=70, c='red', marker='8', alpha=0.5)
ax.scatter(Vmn, Tvmn, s=70, c='red', marker='8', alpha=0.5)
ax.scatter(Vmin, Tvmin, s=70, c='red', marker='8', alpha=0.5)
ax.scatter(Vmax, Tvmax, s=70, c='red', marker='8', alpha=0.5)
ax.text(Vs-5, Tvs+5000, 'Vs')
ax.text(Vmn-5, Tvmn+5000, 'Vmn')
ax.text(Vmin, Tvmin+5000, 'Vmin')
ax.text(Vmax-25, Tvmax+5000, 'Vmax')
ax.set_xlabel('Velocity (ft/s)')
ax.set_ylabel('Power Required (ft·lbs/s)')
ax.grid(True)
plt.show()
