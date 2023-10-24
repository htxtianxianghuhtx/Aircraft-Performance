'''
ISA数据计算器 by HTX
'''

import math

def ISA_calc(h):
    '''
    h 输入的高度 (m)
    
    T0 ISA 海平面温度 (K)
    P0 ISA 海平面气压 (hPa)
    rho0 ISA 海平面大气密度 (kg/m³)

    g ISA 海平面重力加速度 (m/s²)
    M ISA 空气摩尔质量 (kg/mol)
    R 通用气体常数 (J/(mol·K))

    Hlow ISA 对流层-层底高度 (m)
    H1 ISA 对流层-层顶高度/对流层顶-层底高度 (m)
    H2 ISA 对流层顶-层顶高度/第一平流层-层底高度 (m)
    H3 ISA 第一平流层-层顶高度/第二平流层-层底高度 (m)
    H4 ISA 第二平流层-层顶高度/平流层顶-层底高度 (m)
    H5 ISA 平流层顶-层顶高度 (m)
    
    温度梯度 +:递增 -:递减
    lambda1 ISA 对流层(-610~11000m) 温度梯度 (K/m)
    lambda3 ISA 第一平流层(20000~32000m) 温度梯度 (K/m)
    lambda4 ISA 第二平流层(32000~47000m) 温度梯度 (K/m)
    T1to2 ISA 对流层顶(11000~20000m) 温度 (K)
    T3 ISA 第一平流层-层顶(32000m) 温度 (K)
    T4to5 ISA 平流层顶(47000~51000m) 温度 (K)
    P1 ISA 11000m 气压 (hPa)
    P2 ISA 20000m 气压 (hPa)
    P3 ISA 32000m 气压 (hPa)
    P4 ISA 47000m 气压 (hPa)
    rho1 ISA 11000m 大气密度 (kg/m³)
    rho2 ISA 20000m 大气密度 (kg/m³)
    rho3 ISA 32000m 大气密度 (kg/m³)
    rho4 ISA 47000m 大气密度 (kg/m³)
    '''

    T0, P0, rho0 = 288.15, 1013.25, 1.225

    g, M, R = 9.80665, 0.0289644, 8.3144598

    Hlow, H1, H2, H3, H4, H5 = -610, 11000, 20000, 32000, 47000, 51000
    lambda1, lambda3, lambda4 = -0.0065, 0.001, 0.0028
    T1to2 = T0 + lambda1 * H1
    T3 = T1to2 + lambda3 * (H3 - H2)
    T4to5 = T3 + lambda4 * (H4 - H3)
    P1 = P0 * (1 + lambda1 * H1 / T0)**(-g * M / (R * lambda1))
    P2 = P1 * math.exp(-g * M * (H2 - H1) / (R * T1to2))
    P3 = P2 * (1 + lambda3 * (H3 - H2) / T1to2)**(-g * M / (R * lambda3))
    P4 = P3 * (1 + lambda4 * (H4 - H3)/ T3)**(-g * M / (R * lambda4))
    rho1 = rho0 * (1 + lambda1 * H1 / T0)**(-g * M / (R * lambda1) - 1)
    rho2 = rho1 * math.exp(-g * M * (H2 - H1) / (R * T1to2))
    rho3 = rho2 * (1 + lambda3 * (H3 - H2) / T1to2)**(-g * M / (R * lambda3) - 1)
    rho4 = rho3 * (1 + lambda4 * (H4 - H3) / T3)**(-g * M / (R * lambda4) - 1)

    if h < Hlow:
        print('Error: h should not lower than -610m')
        print('T P rho will be -1')
        T = -1
        P = -1
        rho = -1
    elif h <= H1:
        T = T0 + lambda1 * h
        P = P0 * (1 + lambda1 * h / T0)**(-g * M / (R * lambda1))
        rho = rho0 * (1 + lambda1 * h / T0)**(-g * M / (R * lambda1) - 1)
    elif h <= H2:
        T = T1to2
        P = P1 * math.exp(-g * M * (h - H1) / (R * T1to2))
        rho = rho1 * math.exp(-g * M * (h - H1) / (R * T1to2))
    elif h <= H3:
        T = T1to2 + lambda3 * (h - H2)
        P = P2 * (1 + lambda3 * (h - H2) / T1to2)**(-g * M / (R * lambda3))
        rho = rho2 * (1 + lambda3 * (h - H2) / T1to2)**(-g * M / (R * lambda3) - 1)
    elif h <= H4:
        T = T3 + lambda4 * (h - H3)
        P = P3 * (1 + lambda4 * (h - H3)/ T3)**(-g * M / (R * lambda4))
        rho = rho3 * (1 + lambda4 * (h - H3) / T3)**(-g * M / (R * lambda4) - 1)
    elif h <= H5:
        T = T4to5
        P = P4 * math.exp(-g * M * (h - H4) / (R * T4to5))
        rho = rho4 * math.exp(-g * M * (h - H4) / (R * T4to5))
    else:
        print('Error: h should not higher than 51000m')
        print('T P rho will be -2')
        T = -2
        P = -2
        rho = -2
    
    return T, P, rho

print(ISA_calc(11582.4)) #调用函数(11582.4m = 38000ft)
