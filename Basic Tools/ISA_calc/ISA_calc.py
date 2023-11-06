'''
公式
(T P Rho a g)

对于温度恒定区域:
温度 T(h) = T_base (恒定)
气压 P(h) = P_base * e^(-g0 * (h - H_base) / (R_special * T_base))
空气密度 Rho(h) = Rho_base * e^(-g0 * (h - H_base) / (R_special * T_base))

对于温度非恒定区域:
温度 T(h) = T_base + Lambda * (h - H_base)
气压 P(h) = P_base * (1 + Lambda * h / T_base)^(-g0 / (R_special * Lambda))
空气密度 Rho(h) = Rho_base * (1 + Lambda * h / T_base)^(-g0 / (R_special * Lambda) - 1)

a(h) = math.sqrt(Gamma * R_special * T)
g(h) = g0 * (Re / (Re + h))**2

T(h): ISA模型下 给定高度 的 温度 ( 开尔文 | K )
P(h): ISA模型下 给定高度 的 气压 ( 百帕 | hPa )
Rho(h): ISA模型下 给定高度 的 大气密度 ( 千克/立方米 | kg/m³ )
a(h): ISA模型下 给定高度 的 音速 ( 米/秒 | m/s )
g(h): ISA模型下 给定高度 的 重力加速度 ( 米/平方秒 | m/s² )

R: 通用气体常数 ( 焦耳/(摩尔·开尔文) | J/(mol·K) )
M: 空气摩尔质量-干燥空气 ( 千克/摩尔 | kg/mol )
R_special: 特定气体常数-干燥空气 ( 焦耳/(千克·开尔文) | J/(kg·K) )
Re: 地球平均半径 ( 米 | m )
Gamma: 绝热指数-干燥空气 ( 无量纲 )

Lambda: 该段大气层 气温梯度 ( 开尔文/米 | K/m )
T_base: ISA模型下 该段大气层层底 的 温度 ( 开尔文 | K )
P_base: ISA模型下 该段大气层层底 的 气压 ( 百帕 | hPa )
Rho_base: ISA模型下 该段大气层层底 的 大气密度 ( 千克/立方米 | kg/m³ )
g0: ISA海平面 重力加速度 ( 米/平方秒 / m/s² )
'''

import math
import numpy as np

def ISA_calc(h):
    ''' ( 输入 | Input )
    h: 给定高度 ( 米 | m )
    '''
    ''' (输出 | Output)
    T: ISA模型下 给定高度 的 气温 ( 开尔文 | K )
    P: ISA模型下 给定高度 的 气压 ( 百帕 | hPa )
    Rho: ISA模型下 给定高度 的 大气密度 ( 千克/立方米 | kg/m³ )
    a: ISA模型下 给定高度 的 音速 ( 米/秒 | m/s )
    g: ISA模型下 给定高度 的 重力加速度 ( 米/平方秒 | m/s² )
    '''
    ''' ( 说明 | Instruction )
    H_lowest: 对流层-层底位势高度 ( 米 | m )
    H0: 海平面高度 ( 米 | m )
    H1: 对流层-层顶位势高度/对流层顶-层底位势高度 ( 米 | m )
    H2: 对流层顶-层顶位势高度/第一平流层-层底位势高度 ( 米 | m )
    H3: 第一平流层-层顶位势高度/第二平流层-层底位势高度 ( 米 | m )
    H4: 第二平流层-层顶位势高度/平流层顶-层底位势高度 ( 米 | m )
    H5: 平流层顶-层顶位势高度/第一中间层-层底位势高度 ( 米 | m )
    H6: 第一中间层-层顶位势高度/第二中间层-层底位势高度 ( 米 | m )
    H7: 第二中间层-层顶位势高度/中间层顶-层底位势高度 ( 米 | m )

    R: 通用气体常数 ( 焦耳/(摩尔·开尔文) | J/(mol·K) )
    M: 空气摩尔质量-干燥空气 ( 千克/摩尔 | kg/mol )
    R_special: 特定气体常数-干燥空气 ( 焦耳/(千克·开尔文) | J/(kg·K) )
    Re: 地球平均半径 ( 米 | m )
    Gamma: 绝热指数-干燥空气 ( 无量纲 )

    T0: ISA海平面 气温 ( 开尔文 | K )
    P0: ISA海平面 气压 ( 百帕 | hPa )
    Rho0: ISA海平面 大气密度 ( 千克/立方米 / kg/m³ )
    g0: ISA海平面 重力加速度 ( 米/平方秒 / m/s² )

    Lambda1: 对流层 气温梯度 ( 开尔文/米 | K/m )
    Lambda3: 第一平流层 气温梯度 ( 开尔文/米 | K/m )
    Lambda4: 第二平流层 气温梯度 ( 开尔文/米 | K/m )
    Lambda6: 第一中间层 气温梯度 ( 开尔文/米 | K/m )
    Lambda7: 第二中间层 气温梯度 ( 开尔文/米 | K/m )

    T_1to2: 对流层-层顶/对流层顶/第一平流层-层底 气温 ( 开尔文 | K )
    T3: 第一平流层-层顶/第二平流层-层底 气温 ( 开尔文 | K )
    T4to5: 第二平流层-层顶/平流层顶/第一中间层-层底 气温 ( 开尔文 | K )
    T6: 第一中间层-层顶/第二中间层-层底 气温 ( 开尔文 | K )

    P1: 对流层-层顶/对流层顶-层底 气压 ( 百帕 | hPa )
    P2: 对流层顶-层顶/第一平流层-层底 气压 ( 百帕 | hPa )
    P3: 第一平流层-层顶/第二平流层-层底 气压 ( 百帕 | hPa )
    P4: 第二平流层-层底-层顶/平流层顶-层底 气压 ( 百帕 | hPa )
    P5: 平流层顶-层顶/第一中间层-层底 气压 ( 百帕 | hPa )
    P6: 第一中间层-层顶/第二中间层-层底 气压 ( 百帕 | hPa )

    Rho1: 对流层-层顶/对流层顶-层底 气压 ( 千克/立方米 | kg/m³ )
    Rho2: 对流层顶-层顶/第一平流层-层底 气压 ( 千克/立方米 | kg/m³ )
    Rho3: 第一平流层-层顶/第二平流层-层底 气压 ( 千克/立方米 | kg/m³ )
    Rho4: 第二平流层-层底-层顶/平流层顶-层底 气压 ( 千克/立方米 | kg/m³ )
    Rho5: 平流层顶-层顶/第一中间层-层底 气压 ( 千克/立方米 | kg/m³ )
    Rho6: 第一中间层-层顶/第二中间层-层底 气压 ( 千克/立方米 | kg/m³ )
    '''
    
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

print('-610m:', ISA_calc(-610))
print('0m:', ISA_calc(0))
print('11000m:', ISA_calc(11000))
print('11582.4m(38000ft):', ISA_calc(11582.4))
print('20000m:', ISA_calc(20000))
print('32000m:', ISA_calc(32000))
print('47000m:', ISA_calc(47000))
print('51000m:', ISA_calc(51000))
print('71000m:', ISA_calc(71000))
print('84252m:', ISA_calc(84852))
