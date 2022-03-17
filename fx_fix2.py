# 说明：
# 默认atom1为受力原子或交点
# 序号为从0开始的列表号，编号为原子的从1开始的编号
# 论文跑出的工程文件有2个G，如需了解详情联系陈哥 951581823@QQ.com
import math
import os
import random
from copy import deepcopy

import numpy as np
import pandas as pd

# w.write(
#     Atom[0] + '     ' + "{:>2d}".format(int(Atom[1])) + '  ' + "{:<3s}".format(Atom[2]) + ' ' + Atom[2] + ' ' + Atom[
#         4] + ' ' + Atom[4] + '      ' + "{:>.3f}".format(float(Atom[6])) + '  ' + "{:>.3f}".format(
#         float(Atom[7])) + '  ' + "{:>.3f}".format(float(Atom[8])) + '  ' + '1.00' + "{:>.2f}".format(
#         float(Atom[9])) + '           ' + Atom[10])
#
# 全局变量
pi = 3.1415926535898
length = 55 * 2
rp = 8.85418 * 10 ** -12  # relative permittivity 介电常数
e = 1.60217663410 * 10 ** -19
zero = np.array([0.0, 0.0, 0.0]).astype(np.float64)  # 定义零向量
# canshu_vdw = open('F:/毕业设计库/vdW.txt', 'w')
basic_path = 'F:/毕业设计库/basic_data/'  # 基础文件路径

vdw = np.loadtxt('{}vdW.txt'.format(basic_path), dtype=str)
# list = pd.read_csv('F:/毕业设计库/cccc.txt',index_col=0,sep=' ')
ccr = pd.read_csv('{}ccr.txt'.format(basic_path), index_col=0, sep=' ')
rcc = pd.read_csv('{}rcc.txt'.format(basic_path), index_col=0, sep=' ')
# print(vdw)
# print(ccr)
# print(rcc)
# print(ccr.loc['a5', 'D'])
# dih_cccc为四个原子均为α C的二面角参数
dih_cccc = [65.916,  # a1   0
            0.869,  # a2   1
            30.12,  # a3   2
            -0.946,  # a4   3
            -7.939,  # a5   4
            -0.990,  # b1   5
            1.788,  # b2   6
            3.642,  # b3   7
            1.563,  # b4   8
            1.451,  # b5   9
            1.890,  # c1   10
            -0.148,  # c2   11
            0.678,  # c3   12
            0.088,  # c4   13
            -0.516]  # c5   14
# dih_rccc为四个原子中有一个a C的二面角参数
dih_rccc = [[1.522, 1.788],  # a1  0
            [2.155, -94.409],  # a2  1
            [3.077, 94.693],  # a3  2
            [0.000, 86.474],  # a4  3
            [-1.407, -2.573],  # b1  4
            [1.971, 1.363],  # b2  5
            [-1.559, 1.348],  # b3  6
            [0.000, 50.209],  # b4  7
            [1.018, 1.013],  # c1  8
            [1.129, -1.154],  # c2  9
            [7.672, 1.128],  # c3  10
            [1, 28.639]]  # c4  11   # c4原为0，但是无法除以0，所以改为1，反正a=0
#                H_H     H_S    H_T    S_H    S_S     S_T     T_H    T_S    T_T
dih_special = [[-1.938, 3.435, 2.097, 1.274, 1.516, -1.313, -0.257, 0.885, 0.840],  # a1
               [3.945, 2.842, -1.232, 4.331, 4.190, 3.780, 3.855, 2.282, 160.772],  # a2
               [1.592, 2.331, 60.507, 0.745, -1.573, 1.546, -1.361, -1.218, -324.201],  # a3
               [0.607, -1.409, -58.060, -2.017, 1.631, 8.208, 0.904, 3.599, 167.500],  # a4
               [0.866, -1.553, -0.722, -0.903, 2.393, 0.899, 2.033, -0.458, 1.519],  # b1
               [0.477, 2.420, 0.905, 0.287, -0.258, 2.050, 5.918, 3.109, -0.266],  # b2
               [-0.754, -0.723, -0.785, 2.456, 0.823, -0.879, 0.879, 0.845, -0.250],  # b3
               [2.393, 0.839, -0.796, 0.866, -0.843, -9.538, -0.542, -1.498, -0.234],  # b4
               [0.218, 4.282, -0.941, 0.606, 0.965, 0.226, -0.250, 0.876, -0.365],  # c1
               [9.433, 1.467, 0.205, 5.429, 4.273, 3.194, 28.748, 2.538, 2.054],  # c2
               [0.786, 0.938, -3.642, 0.477, 0.214, -0.834, -0.209, -0.187, 2.249],  # c3
               [-0.650, 0.151, 3.398, -0.239, 0.741, 6.577, 0.818, 4.275, 2.494]]  # c4


def angle_emj(atom3, atom1, atom2, atom4):
    # 以AB为交线的 面ABC 和 面ABD的夹角 (C,A,B,D) 对应实参（A，B，C，D)对于形参就是以BC为交线【如需更改，只需交换下形参顺序】
    # 注释方法不可判断补角
    # x1 = float(atom1[6])
    # y1 = float(atom1[7])
    # z1 = float(atom1[8])
    # x2 = float(atom2[6])
    # y2 = float(atom2[7])
    # z2 = float(atom2[8])
    # x3 = float(atom3[6])
    # y3 = float(atom3[7])
    # z3 = float(atom3[8])
    # x4 = float(atom4[6])
    # y4 = float(atom4[7])
    # z4 = float(atom4[8])
    # # a向量
    # ax = x1 - x2
    # ay = y1 - y2
    # az = z1 - z2
    # # b向量
    # bx = x1 - x3
    # by = y1 - y3
    # bz = z1 - z3
    # # c向量
    # cx = x1 - x4
    # cy = y1 - y4
    # cz = z1 - z4
    # # (b,a,c)
    # a = (float(ax) * float(bz) - float(bx) * float(az)) / (float(bx) * float(ay) - float(ax) * float(by))
    # b = -(float(by) * a + float(bz)) / float(bx)
    # c = 1
    # # (e,d,f)
    # # (e,d,f)
    # d = (float(ax) * float(cz) - float(cx) * float(az)) / (float(cx) * float(ay) - float(ax) * float(cy))
    # e = -(float(cy) * float(d) + float(cz)) / float(cx)
    # f = 1
    # g = (b ** 2 + a ** 2 + c ** 2) ** 0.5
    # h = (e ** 2 + d ** 2 + f ** 2) ** 0.5
    # cosC = (b * e + a * d + c * f) / (g * h)

    # 向量AB
    AB = np.array(
        [float(atom2[2]) - float(atom1[2]), float(atom2[3]) - float(atom1[3]), float(atom2[4]) - float(atom1[4])])
    # 向量AC
    AC = np.array(
        [float(atom3[2]) - float(atom1[2]), float(atom3[3]) - float(atom1[3]), float(atom3[4]) - float(atom1[4])])
    # 向量AD
    AD = np.array(
        [float(atom4[2]) - float(atom1[2]), float(atom4[3]) - float(atom1[3]), float(atom4[4]) - float(atom1[4])])
    # # 向量BD
    # BD = np.array(
    #     [float(atom4[2]) - float(atom2[2]), float(atom4[3]) - float(atom2[3]), float(atom4[4]) - float(atom2[4])])
    # # 向量BC
    # BC = np.array(
    #     [float(atom3[2]) - float(atom2[2]), float(atom3[3]) - float(atom2[3]), float(atom3[4]) - float(atom2[4])])

    # 向量ED——D垂AB于点E
    ED = AD - AB * (np.linalg.norm(AB) * math.cos(angle_jj(atom1, atom2, atom4) / np.linalg.norm(AB)))
    # 向量FC——C垂AB于点F
    FC = AC - AB * (np.linalg.norm(AC) * math.cos(angle_jj(atom1, atom2, atom3) / np.linalg.norm(AB)))

    # 向量ED和向量FC的夹角即为二面角
    if np.linalg.norm(ED) * np.linalg.norm(FC) == 0:
        cosC = 1
    else:
        cosC = np.dot(ED, FC) / (np.linalg.norm(ED) * np.linalg.norm(FC))

    angle = math.acos(float(cosC))
    return angle


def angle_jj(atom1, atom2, atom3):  # 以Atom1为交点

    # x1 = float(atom1[6])
    # y1 = float(atom1[7])
    # z1 = float(atom1[8])
    # x2 = float(atom2[6])
    # y2 = float(atom2[7])
    # z2 = float(atom2[8])
    # x3 = float(atom3[6])
    # y3 = float(atom3[7])
    # z3 = float(atom3[8])
    #
    # a1 = x2 - x1
    # a2 = y2 - y1
    # a3 = z2 - z1
    # b1 = x3 - x1
    # b2 = y3 - y1
    # b3 = z3 - z1

    # a = np.array([a1, a2, a3])
    # b = np.array([[b1], [b2], [b3]])

    a = np.array(
        [float(atom2[2]) - float(atom1[2]), float(atom2[3]) - float(atom1[3]), float(atom2[4]) - float(atom1[4])])
    b = np.array(
        [[float(atom3[2]) - float(atom1[2])], [float(atom3[3]) - float(atom1[3])], [float(atom3[4]) - float(atom1[4])]])

    cosC = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

    # a = ((float(x1) - float(x2)) ** 2 + (float(y1) - float(y2)) ** 2 + (float(z1) - float(z2)) ** 2) ** 0.5
    # b = ((float(x2) - float(x3)) ** 2 + (float(y2) - float(y3)) ** 2 + (float(z2) - float(z3)) ** 2) ** 0.5
    # c = ((float(x1) - float(x3)) ** 2 + (float(y1) - float(y3)) ** 2 + (float(z1) - float(z3)) ** 2) ** 0.5
    # cosC = (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)

    angle = math.acos(cosC)
    return angle


# 以下三个函数为判断键长函数
# def c(atom):  # 判断是否为 α C，是则返回1 否则返回0
#     a = "_C" in atom[1]
#     if a == 'true':
#         b = 1
#     else:
#         b = 0
#     return b


# def rbound(atom):
#     if '_R' in atom[1]:
#         return atom[11]
#     else:
#         return 'false'


def lbound(atom1, atom2):
    a = float(atom1[6])
    b = float(atom2[6])
    if a == b:
        return a
    elif a == 3.8:
        return b
    elif b == 3.8:
        return a
    else:
        return False


def Normal_vector(atom1, atom2, atom3):  # 计算法向量 不确定方向型 若三点共线，则返回一个垂直于改线的法向量
    # 方法，解齐次线性方程组

    x1 = float(atom2[2]) - float(atom3[2])
    y1 = float(atom2[3]) - float(atom3[3])
    z1 = float(atom2[4]) - float(atom1[4])

    x2 = float(atom1[2]) - float(atom3[2])
    y2 = float(atom1[3]) - float(atom3[3])
    z2 = float(atom1[4]) - float(atom3[4])

    # x3 = float(atom1[2]) - float(atom2[2])
    # y3 = float(atom1[3]) - float(atom2[3])
    # z3 = float(atom1[4]) - float(atom2[4])
    #
    # A = np.array([[x1,y1,z1],
    #         [x2,y2,z2],
    #         [x3,y3,z3]])
    # B = np.array([[0],[0],[0]])
    # x=np.linalg.solve(A, B)
    #
    #
    # if x1 != 0 and y2 * x1 - x2 * y1 != 0:
    #     b = (x2 * z1 - z2 * x1) / (y2 * x1 - x2 * y1)
    #     a = -(b * y1 + z1) / x1
    #     c = 1
    # elif y2 * x1 - x2 * y1 == 0:
    #     c = 0
    #     a = y1 / x1
    #     b = -1
    # elif x1 == 0 and y1 == 0:
    #     a = 1
    a = y1 * z2 - y2 * z1
    b = z1 * x2 - z2 * x1
    c = x1 * y2 - x2 * y1

    fx = np.array([a, b, c]).astype(np.float64)

    if fx[0] + fx[1] + fx[2] < 10 ** 5:  # 三个点在一条直线时，力的方向为0向量 设力的方向为垂直于该向量的任意单位向量
        zero_forbid = deepcopy(zero)
        if x1 == 0: zero_forbid[0] = random.randint(1, 100)
        if y1 == 0: zero_forbid[1] = random.randint(1, 100)
        if z1 == 0: zero_forbid[2] = random.randint(1, 100)
        if np.linalg.norm(zero_forbid) == 0:
            zero_forbid[1] = random.randint(1, 100)
            zero_forbid[2] = random.randint(1, 100)
            zero_forbid[0] = -(zero_forbid[1] * y1 + zero_forbid[2] * z1) / x1

        fx = zero_forbid / np.linalg.norm(zero_forbid)
    fx = fx / np.linalg.norm(fx)

    return fx


def mkdir(path):
    folder = os.path.exists(path)

    if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径
        print("---  新目录地址...  ---")
        print("---  已创建新目录  ---")

    else:
        print("---  地址已存在，不再创建  ---")


def distance(atom1, atom2):
    dis = ((float(atom1[2]) - float(atom2[2])) ** 2 + (float(atom1[3]) - float(atom2[3])) ** 2 + (
            float(atom1[4]) - float(atom2[4])) ** 2) ** 0.5
    return dis


def distance_e(atom1, atom2):
    dis = ((float(atom1[9]) - float(atom2[9])) ** 2 + (float(atom1[10]) - float(atom2[10])) ** 2 + (
            float(atom1[11]) - float(atom2[11])) ** 2) ** 0.5
    return dis


def rad():  # 随机生成-1或1
    rad = random.randint(0, 1)
    if rad == 0:
        rad = -1
    return rad


def group(atom):
    if atom[1][0] == 'E' or atom[1][0] == 'A' or atom[1][0] == 'L' or atom[1][0] == 'M' or atom[1][0] == 'Q' or atom[1][
        0] == 'K' or atom[1][0] == 'R' or atom[1][0] == 'H':
        return 'H'
    elif atom[1][0] == 'V' or atom[1][0] == 'I' or atom[1][0] == 'Y' or atom[1][0] == 'C' or atom[1][0] == 'W' or \
            atom[1][0] == 'F' or atom[1][0] == 'T':
        return 'S'
    elif atom[1][0] == 'G' or atom[1][0] == 'N' or atom[1][0] == 'P' or atom[1][0] == 'S' or atom[1][0] == 'D':
        return 'T'
    else:
        return 'error'


# 以下为能量项和力项

def U_bound(atom1, atom2):  # 设定键能常量为325   1250过大  100的过小
    r = distance(atom1, atom2)
    U = 325 * (r - lbound(atom1, atom2)) ** 2
    return U


def force_bound(atom1, atom2):  # 键力矢量 作用于atom1
    l = lbound(atom1, atom2)
    r = distance(atom1, atom2)
    f1 = np.array(
        [float(atom1[2]) - float(atom2[2]), float(atom1[3]) - float(atom2[3]), float(atom1[4]) - float(atom2[4])])
    f2 = f1 / np.linalg.norm(f1)  # 单位向量指向atom1
    fx = 650 * (l - r) * f2
    fx = fx.astype('float64')
    return fx


def U_vdW(atom1, atom2):  # 主链碳原子均被粗粒化为Glycine原子
    c1 = c2 = d1 = d2 = 0
    r = distance(atom1, atom2)
    for i in range(len(vdw)):
        if vdw[i][0] in atom1[1]:
            c1 = float(vdw[i][1])
            d1 = float(vdw[i][2])
            break
    for i in range(len(vdw)):
        if vdw[i][0] in atom2[1]:
            c2 = float(vdw[i][1])
            d2 = float(vdw[i][2])
            break

    c = (c1 + c2) / 2
    d = math.sqrt(d1 * d2)
    U = 4 * 4.18585 * d * ((c / r) ** 12 - (c / r) ** 6)

    return U


def force_vdW(atom1, atom2):  # 范德华力矢量  指向atom1
    c1 = c2 = d1 = d2 = 0
    for i in range(len(vdw)):
        if vdw[i][0] in atom1[1]:
            c1 = float(vdw[i][1])
            d1 = float(vdw[i][2])
            break
    for i in range(len(vdw)):
        if vdw[i][0] in atom2[1]:
            c2 = float(vdw[i][1])
            d2 = float(vdw[i][2])
            break

    r = distance(atom1, atom2)
    c = (c1 + c2) / 2
    d = math.sqrt(d1 * d2)
    f = 4 * 4.18585 * d * (-12 * (c / r) ** 12 + 6 * (c / r) ** 6) / r ** 2
    # 力的方向，指向atom1
    f1 = np.array(
        [float(atom1[2]) - float(atom2[2]), float(atom1[3]) - float(atom2[3]), float(atom1[4]) - float(atom2[4])])

    f2 = - f1 / np.linalg.norm(f1)
    fx = f * f2

    return fx


# def U_angle(atom1, atom2, atom3):  # 键角能 A为交点 B为作用点  只做用于R基，主函数中做判断  计算两个能：C-C-R 和 R-C-C，C-C-C不考虑
#     x = angle_jj(atom1, atom2, atom3)
#     # U = 3 * math.exp(-((θ - 3) / 3) ** 2)
#     #U = 0.5 * 25 * (math.cos(θ) + 1)
#     record = int(atom2[0]) - 1  #记录作用点序号
#     if '_R' in atom1[1]:
#         ATOM[record -1],ATOM[record],ATOM[record -3]
#         ATOM[record - 1], ATOM[record], ATOM[record + 1]
#     while (i < 5):
#         U += dih_cccc[i] * math.exp(-((x - dih_cccc[i + 5]) / dih_cccc[i + 10]) ** 2)
#         i += 1
#     return U
# print(ccr['D'][2])
def U_angle(atom1, atom2, atom3, atom4):  # 键角能  只做用于R基(atom3)，主函数中做判断  计算两个能：C-C-R 和 R-C-C，C-C-C不考虑

    U = 0
    x1 = 0
    x2 = 0
    record = int(atom3[12]) - 1  # 记录作用点序号
    acid = atom3[1][0]  # 记录氨基酸单字代号
    # ccr_g = ccr.loc[acid]
    # rcc_g = rcc.loc[acid]

    if '_R' in atom3[1]:
        if record != 1:  # 第一个R级无C-C-R
            x1 = angle_jj(atom2, atom3, atom1)  # C-C-R
        if record != length - 1:  # 最后一个R基无R-C-C
            x2 = angle_jj(atom2, atom3, atom4)  # R-C-C
        for i in range(5):
            if x1 != 0:
                U += ccr[acid][i] * math.exp(-(float(x1 - ccr[acid][i + 5]) / float(ccr[acid][i + 10])) ** 2)  # C-C-R
            if x2 != 0:
                # print(rcc[acid][i + 10])
                U += rcc[acid][i] * math.exp(-(float(x2 - rcc[acid][i + 5]) / float(rcc[acid][i + 10])) ** 2)  # R-C-C
                # print(float(rcc[acid][i + 10]))
        return U * 4.1858518
    else:
        return 'not R'


def U_angle_fix(atom1, atom2, atom3, atom4, atom5):
    x1 = angle_jj(atom1, atom2, atom3)  # rcc角
    x2 = angle_jj(atom4, atom5, atom3)  # ccr角
    U = 0
    record = int(atom3[12]) - 1  # 记录计算点的序号
    acid_rcc = atom2[1][0]
    acid_ccr = atom5[1][0]
    for i in range(5):
        if record != 108 and acid_ccr != 'G':  # 最后一个原子无CCR
            U += ccr[acid_ccr][i] * math.exp(
                -(float(x1 - ccr[acid_ccr][i + 5]) / float(ccr[acid_ccr][i + 10])) ** 2)  # C-C-R
        if record != 0 and acid_rcc != 'G':  # 第一个一个原子无RCC
            U += rcc[acid_rcc][i] * math.exp(
                -(float(x2 - rcc[acid_rcc][i + 5]) / float(rcc[acid_rcc][i + 10])) ** 2)  # R-C-C
    return U * 4.1858518


def force_angle_fix(atom1, atom2, atom3, atom4, atom5):
    x1 = angle_jj(atom1, atom2, atom3)  # rcc角
    x2 = angle_jj(atom4, atom5, atom3)  # ccr角
    acc_ccr = acc_rcc = 0
    record = int(atom3[12]) - 1  # 记录计算点的序号
    acid_rcc = atom2[1][0]
    acid_ccr = atom5[1][0]
    for i in range(5):
        if record != 108 and acid_ccr != 'G':
            acc_ccr += - 2 * (x1 - ccr[acid_ccr][i + 5]) * (  # pandas 调用是test[列][行]
                    ccr[acid_ccr][i] * math.exp(-(float(x1 - ccr[acid_ccr][i + 5]) / ccr[acid_ccr][i + 10]) ** 2)) / \
                       ccr[acid_ccr][
                           i + 10] ** 2  # C-C-R
        if record != 0 and acid_rcc != 'G':
            acc_rcc += - 2 * (x2 - rcc[acid_rcc][i + 5]) * (
                    rcc[acid_rcc][i] * math.exp(-(float(x1 - rcc[acid_rcc][i + 5]) / rcc[acid_rcc][i + 10]) ** 2)) / \
                       rcc[acid_rcc][
                           i + 10] ** 2  # R-C-C
    # 指向前一个R   R-C-C
    f1 = np.array(  #
        [float(atom2[2]) - float(atom3[2]), float(atom2[3]) - float(atom3[3]),
         float(atom2[4]) - float(atom3[4])]).astype(np.float64)
    f_rcc = f1 / np.linalg.norm(f1)
    # 指向后一个R   C-C-R
    f2 = np.array(
        [float(atom5[2]) - float(atom3[2]), float(atom5[3]) - float(atom3[3]),
         float(atom5[4]) - float(atom3[4])]).astype(np.float64)
    f_ccr = f2 / np.linalg.norm(f2)

    if f_ccr[0] + f_ccr[0] + f_ccr[0] < 10 ** 5:
        f_ccr = Normal_vector(atom1, atom2, atom3)
    if f_rcc[0] + f_rcc[0] + f_rcc[0] < 10 ** 5:
        f_rcc = Normal_vector(atom3, atom4, atom5)

    f = acc_rcc * f_rcc + acc_ccr * f_ccr
    return f * 4.1858518


# def force_angle(atom1, atom2, atom3):  # 虚拟键角能  交点为atom1，作用于atom2，方向C > B 作用强度为9
#     x = angle_jj(atom1, atom2, atom3)
#     f = -2 * (x - 3) * math.exp(-((x - 3) / 3) ** 2)
#     dir = np.array(
#         [float(atom3[2]) - float(atom2[2]), float(atom3[3]) - float(atom2[3]),
#          float(atom3[4]) - float(atom2[4])]).astype(np.float64)
#     fx = dir / np.linalg.norm(dir)
#     if f == 0:
#         return 0
#     elif f < 0:
#         return dir * f
#     else:
#         return -1 * dir * f

def force_angle(atom1, atom2, atom3, atom4):  # 虚拟键角能  交点为atom1，作用于atom2，方向C > B
    acc_rcc = 0
    acc_ccr = 0
    x1 = 0
    x2 = 0
    record = int(atom3[12]) - 1  # 记录作用点序号
    acid = atom3[1][0]  # 记录氨基酸单字代号
    # ccr_g = ccr.loc[acid]
    # rcc_g = rcc.loc[acid]

    if '_R' in atom3[1]:
        if record != 1:  # 第一个R级无C-C-R
            x1 = angle_jj(atom2, atom3, atom1)  # C-C-R
        if record != length - 1:  # 最后一个R基无R-C-C
            x2 = angle_jj(atom2, atom3, atom4)  # R-C-C
        for i in range(5):
            if x1 != 0:
                acc_ccr += - 2 * (x1 - ccr[acid][i + 5]) * (  # pandas 调用是test[列][行]
                        ccr[acid][i] * math.exp(-(float(x1 - ccr[acid][i + 5]) / ccr[acid][i + 10]) ** 2)) / ccr[acid][
                               i + 10] ** 2  # C-C-R
            if x2 != 0:
                acc_rcc += - 2 * (x2 - rcc[acid][i + 5]) * (
                        rcc[acid][i] * math.exp(-(float(x1 - rcc[acid][i + 5]) / rcc[acid][i + 10]) ** 2)) / rcc[acid][
                               i + 10] ** 2  # R-C-C
    else:
        print('not R')
    # 指向前一个C   C-C-R
    f1 = np.array(  #
        [float(atom1[2]) - float(atom3[2]), float(atom1[3]) - float(atom3[3]),
         float(atom1[4]) - float(atom3[4])]).astype(np.float64)
    f_ccr = f1 / np.linalg.norm(f1)
    # 指向后一个C   R-C-C
    f2 = np.array(
        [float(atom4[2]) - float(atom3[2]), float(atom4[3]) - float(atom3[3]),
         float(atom4[4]) - float(atom3[4])]).astype(np.float64)
    f_rcc = f2 / np.linalg.norm(f2)

    if f_ccr[0] + f_ccr[0] + f_ccr[0] < 10 ** 5:
        f_ccr = Normal_vector(atom2, atom3, atom4)
    if f_rcc[0] + f_rcc[0] + f_rcc[0] < 10 ** 5:
        f_rcc = Normal_vector(atom1, atom2, atom3)

    f = acc_rcc * f_rcc + acc_ccr * f_ccr
    return f * 4.1858518


def U_dih(atom1, atom2, atom3, atom4):  # 以BC为交线的 (A,B,C,D)面ABC 和 面BCD的夹角
    x = angle_emj(atom1, atom2, atom3, atom4)
    if '_R' in atom1[1]:
        count = 1
    else:
        count = 0
    # if '_r' in atom1[1]: count += 1  # 判断四个原子中有多少个R基
    # if '_r' in atom2[1]: count += 1
    # if '_r' in atom3[1]: count += 1
    # if '_r' in atom4[1]: count += 1

    # count_1 = 0
    # count_2 = 0
    # count_3 = 0
    # count_4 = 0
    # count_5 = 0
    # count_6 = 0
    # if 'H_' in atom2[1]: count_1 += 1  # 判断是否中间两个原子同时含有H、S、T
    # if 'S_' in atom2[1]: count_2 += 1
    # if 'T_' in atom2[1]: count_3 += 1
    # if 'H_' in atom3[1]: count_4 += 1
    # if 'S_' in atom3[1]: count_5 += 1
    # if 'T_' in atom3[1]: count_6 += 1
    # count_0 = count_1 + count_2 + count_3 + count_4 + count_5 + count_6

    count_1 = []
    count_2 = []
    title = [[0, 1, 2],
             [3, 4, 5],
             [6, 7, 8]]
    if group(atom2) == 'H': count_1.extend([0])  # 判断是否中间两个原子同时含有H、S、T
    if group(atom2) == 'S': count_1.extend([1])
    if group(atom2) == 'T': count_1.extend([2])
    if group(atom3) == 'H': count_2.extend([0])
    if group(atom3) == 'S': count_2.extend([1])
    if group(atom3) == 'T': count_2.extend([2])
    count_0 = len(count_1) + len(count_2)

    U = 0
    if count == 0:  # R-R-R-R情况
        if count_0 < 2:
            for i in range(5):
                U += dih_cccc[i] * math.exp(-((x - dih_cccc[i + 5]) / dih_cccc[i + 10]) ** 2)
        elif count_0 == 2:  # R-Ri-Rj-R情况
            read = title[count_1[0]][count_2[0]]
            for i in range(4):
                U += dih_special[i][read] * math.exp(-((x - dih_special[i + 4][read]) / dih_special[i + 8][read]) ** 2)

    elif count == 1:  # or count == 2:  # 注释为R-C-C-R不使用
        for i in range(4):  # R-C-C-C
            if 'P_' in atom1[1]:
                U += dih_rccc[i][1] * math.exp(-((x - dih_rccc[i + 4][1]) / dih_rccc[i + 8][1]) ** 2)
            elif 'G_' in atom1[1]:
                U = 0
            else:
                U += dih_rccc[i][0] * math.exp(-((x - dih_rccc[i + 4][0]) / dih_rccc[i + 8][0]) ** 2)
        # if count_0 == 2:  # 对于中间两个碳原子为H、S、T的 采用cccc-special参数
        #     read = title[count_1[0]][count_2[0]]
        #     for i in range(5):
        #         U += dih_special[read][i] * math.exp(-((x - dih_special[read][i + 4]) / dih_special[read][i + 8]) ** 2)



    else:
        U = False
    return U * 4.1858518


def force_dih(atom1, atom2, atom3, atom4):  # 二面角能力 作用于atom4

    # # 作用方向 D > A 会造成R基偏到主链上
    # # dir = np.array([float(atom4[2]) - float(atom1[2]), float(atom4[3]) - float(atom1[3]),
    # #                 float(atom4[4]) - float(atom1[4])])  # 二面角能力方向CD
    # # CB = np.array(
    # #     [float(atom2[2]) - float(atom3[2]), float(atom2[3]) - float(atom3[3]), float(atom2[4]) - float(atom1[4])])
    # x1 = float(atom2[2]) - float(atom3[2])
    # y1 = float(atom2[3]) - float(atom3[3])
    # z1 = float(atom2[4]) - float(atom3[4])
    #
    # # CD = np.array(
    # #     [float(atom4[2]) - float(atom3[2]), float(atom4[3]) - float(atom3[3]), float(atom4[4]) - float(atom3[4])])
    # x2 = float(atom4[2]) - float(atom3[2])
    # y2 = float(atom4[3]) - float(atom3[3])
    # z2 = float(atom4[4]) - float(atom3[4])

    # b = (x2 * z1 - z2 * x1) / (y2 * x1 - x2 * y1)
    # a = -(b * y1 + z1) / x1
    # c = 1
    # f1 = np.array([a,b,c])
    # f2 = f1 / np.linalg.norm(f1)

    f2 = Normal_vector(atom2, atom3, atom4)

    place = deepcopy(atom4)  # 判断方向
    place[2] = place[2] + f2[0] * 0.01
    place[3] = place[3] + f2[1] * 0.01
    place[4] = place[3] + f2[2] * 0.01

    dir = f2.astype(np.float64)
    # if dir[0] + dir[1] + dir[2] < 10 ** 5:  # 三个点在一条直线时，力的方向为0向量 设力的方向为垂直于改向量的任意向量
    #     zero_forbid = np.array([-2 / x1, 1, 1]).astype(np.float64)
    #     dir = zero_forbid / np.linalg.norm(zero_forbid)

    if distance(place, atom1) > distance(atom4, atom1):
        dir = - f2.astype(np.float64)

    x = angle_emj(atom1, atom2, atom3, atom4)
    if '_R' in atom1[1]:
        count = 1
    else:
        count = 0
    # if '_r' in atom1[1]: count += 1  # 判断四个原子中有多少个R基
    # if '_r' in atom2[1]: count += 1
    # if '_r' in atom3[1]: count += 1
    # if '_r' in atom4[1]: count += 1

    # count_1 = 0
    # count_2 = 0
    # count_3 = 0
    # count_4 = 0
    # count_5 = 0
    # count_6 = 0
    # if 'H_' in atom2[1]: count_1 += 1  # 判断是否中间两个原子同时含有H、S、T
    # if 'S_' in atom2[1]: count_2 += 1
    # if 'T_' in atom2[1]: count_3 += 1
    # if 'H_' in atom3[1]: count_4 += 1
    # if 'S_' in atom3[1]: count_5 += 1
    # if 'T_' in atom3[1]: count_6 += 1
    # count_0 = count_1 + count_2 + count_3 + count_4 + count_5 + count_6

    count_1 = []
    count_2 = []
    title = [[0, 1, 2],
             [3, 4, 5],
             [6, 7, 8]]
    if group(atom2) == 'H': count_1.extend([0])  # 判断是否中间两个原子同时含有H、S、T
    if group(atom2) == 'S': count_1.extend([1])  # 可以使用count_1.append(1) 可能效率会更高但影响不大
    if group(atom2) == 'T': count_1.extend([2])
    if group(atom3) == 'H': count_2.extend([0])
    if group(atom3) == 'S': count_2.extend([1])
    if group(atom3) == 'T': count_2.extend([2])
    count_0 = len(count_1) + len(count_2)
    if count_0 == 2:
        read = title[count_1[0]][count_2[0]]
    i = 0
    force = 0

    U = 0
    if count == 0:  # C-C-C-C情况
        if count_0 < 2:  # 以下部分可以注释掉无需判断，但是不影响结果输出。即有冗余代码count_0是一定等于2的
            for i in range(5):
                force += (- 2 * dih_cccc[i] * (x - dih_cccc[i + 5]) * math.exp(
                    -((x - dih_cccc[i + 5]) / dih_cccc[i + 10]) ** 2)) / (dih_cccc[i + 10]) ** 2
        elif count_0 == 2:  # C-Ci-Cj-C 情况
            for i in range(4):
                force += (- 2 * dih_special[i][read] * (x - dih_special[i + 4][read]) * math.exp(
                    -((x - dih_special[i + 4][read]) / dih_special[i + 8][read]) ** 2)) / (
                             dih_special[i + 8][read]) ** 2

    elif count == 1:  # or count == 2:  # R-C-C-C  注释为R-Ci-Cj-C不使用
        for i in range(4):
            if 'P_' in atom1[1]:
                force += (- 2 * dih_rccc[i][1] * (x - dih_rccc[i + 4][1]) * math.exp(
                    -((x - dih_rccc[i + 4][1]) / dih_rccc[i + 8][1]) ** 2)) / (dih_rccc[i + 8][1]) ** 2
            elif 'G_' in atom1[1]:
                force = 0
            else:
                force += (- 2 * dih_rccc[i][0] * (x - dih_rccc[i + 4][0]) * math.exp(
                    -((x - dih_rccc[i + 4][0]) / dih_rccc[i + 8][0]) ** 2)) / (dih_rccc[i + 8][0]) ** 2

        # else:
        #     for i in range(4):
        #         force += (- 2 * dih_special[i][read] * (x - dih_special[i + 4][read]) * math.exp(
        #             -((x - dih_special[i + 4][read]) / dih_special[i + 8][read]) ** 2)) / (
        #                      dih_special[i + 8][read]) ** 2


    # if count == 0:
    #     if count_0 < 2:
    #         while i < 5:
    #             force += (- 2 * dih_cccc[i] * (x - dih_cccc[i + 5]) * math.exp(
    #                 -((x - dih_cccc[i + 5]) / dih_cccc[i + 10]) ** 2)) / (dih_cccc[i + 10]) ** 2
    #             i += 1
    #     elif count_2 == 2:
    #         read = title[count_1[0]][count_2[0]]
    #         while i < 4:
    #             force += (- 2 * dih_special[i][read] * (x - dih_special[i + 4][read]) * math.exp(
    #                 -((x - dih_special[i + 4][read]) / dih_special[i + 8][read]) ** 2)) / (
    #                          dih_special[i + 8][read]) ** 2
    #             i += 1
    #
    # elif count == 1 or count == 2:  # 由于原文未给出有两个R基的情况，故套用一个R的情况
    #     if count_2 == 2:  # 同理对于中间两个碳原子为H、S、T的 采用cccc-special参数
    #         read = title[count_1[0]][count_2[0]]
    #         while i < 4:
    #             force += (- 2 * dih_special[i][read] * (x - dih_special[i + 4][read]) * math.exp(
    #                 -((x - dih_special[i + 4][read]) / dih_special[i + 8][read]) ** 2)) / (
    #                          dih_special[i + 8][read]) ** 2
    #             i += 1
    #     else:
    #         while i < 4:
    #             force += (- 2 * dih_cccc[i] * (x - dih_cccc[i + 5]) * math.exp(
    #                 -((x - dih_cccc[i + 5]) / dih_cccc[i + 10]) ** 2)) / (dih_cccc[i + 10]) ** 2
    #             i += 1
    else:
        force = False

    fx = force * dir
    return fx * 4.1858518


def U_ele(atom1, atom2):  # 电荷能
    U = e ** 2 * 6.02 * 10 ** 23 * 10 ** -3 * float(atom1[7]) * float(atom2[7]) / (
            4 * pi * rp * 20 * 10 ** -10 * distance_e(atom1, atom2))
    return U


def force_ele(atom1, atom2):  # 电荷力矢量 B指向A 作用于atom1
    F = U_ele(atom1, atom2) / (distance_e(atom1, atom2))
    # 电荷力方向
    f1 = np.array(
        [float(atom1[2]) - float(atom2[2]), float(atom1[3]) - float(atom2[3]), float(atom1[4]) - float(atom2[4])])
    f2 = f1 / np.linalg.norm(f1)
    fx = f2 * F
    return fx
