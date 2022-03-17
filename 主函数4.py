import time

from fx_fix2 import *

# 蛋白质 4DCZ
# 以下为RCSP上的序列 共计92个氨基酸
# protein = ['K_', 'Q_', 'E_', 'Q_', 'P_', 'E_', 'I_', 'N_', 'L_', 'D_', 'H_', 'V_', 'V_', 'E_', 'Q_', 'T_', 'I_', 'K_',
#            'K_', 'V_', 'Q_', 'Q_', 'N_', 'Q_', 'N_', 'Q_', 'N_', 'K_', 'D_', 'P_', 'D_', 'E_', 'L_', 'R_', 'S_', 'K_',
#            'V_', 'P_', 'G_', 'E_', 'V_', 'T_', 'A_', 'S_', 'D_', 'W_', 'E_', 'A_', 'L_', 'V_', 'G_', 'D_', 'T_', 'R_',
#            'Y_', 'G_', 'Y_', 'F_', 'D_', 'E_', 'T_', 'G_', 'D_', 'W_', 'S_', 'W_', 'K_', 'G_', 'Y_', 'F_', 'D_', 'E_',
#            'Q_', 'G_', 'K_', 'W_', 'V_', 'W_', 'N_', 'E_', 'P_', 'V_', 'D_', 'S_', 'L_', 'E_', 'H_', 'H_', 'H_', 'H_',
#            'H_', 'H_']

# 以下为pdb中的蛋白序列 共计55个氨基酸 为26-80的切片
# protein = protein[25:79]  # 具体如下
protein = ['Q_', 'N_', 'K_', 'D_', 'P_', 'D_', 'E_', 'L_', 'R_', 'S_', 'K_', 'V_', 'P_', 'G_', 'E_', 'V_', 'T_', 'A_',
           'S_', 'D_', 'W_', 'E_', 'A_', 'L_', 'V_', 'G_', 'D_', 'T_', 'R_', 'Y_', 'G_', 'Y_', 'F_', 'D_', 'E_', 'T_',
           'G_', 'D_', 'W_', 'S_', 'W_', 'K_', 'G_', 'Y_', 'F_', 'D_', 'E_', 'Q_', 'G_', 'K_', 'W_', 'V_', 'W_', 'N_',
           'E_']

localtime = time.strftime("%Y_%m_%d", time.localtime())  # "%a %b %d %H:%M:%S %Y"完整格式，如需准确时间可拓展
m_s = time.strftime("%H_%M", time.localtime())  # "%M:%S"小时和分钟，区分每次测试
time_all = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())  # 全部时间

path = 'F:/毕业设计库/测试数据/{:s}/{}/4/'.format(localtime, m_s)  # 路径地址，自定义
mkdir(path)  # 创建路径
# 读取20个氨基酸粗粒文件
pl = np.loadtxt('F:/毕业设计库/protein_ele.txt', dtype=str)  # protein list
rmsd_norm = np.loadtxt('F:/毕业设计库/标准蛋白.txt', dtype=str)  # 标准蛋白
culi = open('{}粗粒化后.txt'.format(path), 'w')
pdb = open('{}粗粒化后.pdb'.format(path), 'w')
CON = open('{}CONECT.txt'.format(path), 'w')
sy = open('{}说明.txt'.format(path), 'w')
note = open('{}note.txt'.format(path), 'w')
note.write('tome' + ' ' + 'step' + ' ' + 'm' + ' ' + 'energy' + ' ' + 'RMSD' + '\n')
sy.write(time_all + '\n' + 'Description :' + '\n')
sy.close()
# result = open('F:/毕业设计库/result.txt', 'w')
# pl = pd.read_csv('F:/毕业设计库/粗粒数据.txt', sep=' ')
i = 0
k = 1

ATOM = []

for i in range(len(protein)):
    for j in range(len(pl)):
        if '_C' in pl[j][1] and protein[i] in pl[j][1]:
            z = (-1) ** 0.5
            while isinstance(z, complex):  # 判断z是否为复数
                x = random.uniform(0, 3.8)
                y = random.uniform(-1.8, 1.8)
                z = (3.8 ** 2 - x ** 2 - y ** 2) ** 0.5

            if i == 0:
                ATOM.extend([
                    [pl[j][0], pl[j][1], 0, 0, 0, pl[j][2], float(pl[j][3]),
                     pl[j][4],
                     float(pl[j][5]),
                     0 + rad() * float(pl[j][5]) ** 0.333, 0 + rad() * float(pl[j][5]) ** 0.333,
                     0 + rad() * float(pl[j][5]) ** 0.333, k]])
            else:
                ATOM.extend([
                    [pl[j][0], pl[j][1], ATOM[k - 3][2] + x, ATOM[k - 3][3] + y, ATOM[k - 3][4] + z, pl[j][2],
                     float(pl[j][3]),
                     pl[j][4],
                     float(pl[j][5]),
                     ATOM[k - 3][2] + x + rad() * float(pl[j][5]) ** 0.333,
                     ATOM[k - 3][3] + y + rad() * float(pl[j][5]) ** 0.333,
                     ATOM[k - 3][4] + z + rad() * float(pl[j][5]) ** 0.333, k]])

            ATOM.extend([
                [pl[j + 1][0], pl[j + 1][1], ATOM[k - 1][2] + rad() * float(pl[j + 1][3]) ** 0.333,
                 ATOM[k - 1][3] + rad() * float(pl[j + 1][3]) ** 0.333,
                 ATOM[k - 1][4] + rad() * float(pl[j + 1][3]) ** 0.333, pl[j + 1][2], float(pl[j + 1][3]),
                 pl[j + 1][4],
                 float(pl[j + 1][5]),
                 ATOM[k - 1][2] + float(pl[j + 1][3]) + rad() * float(pl[j + 1][5]) ** 0.333,
                 ATOM[k - 1][3] + float(pl[j + 1][3]) + rad() * float(pl[j + 1][5]) ** 0.333,
                 ATOM[k - 1][4] + float(pl[j + 1][3]) + rad() * float(pl[j + 1][5]) ** 0.333, k + 1]])
            k += 2
# 以下为错误生成
# for i in range(len(protein)):
#     j = 0
#     for j in range(len(pl)):
#         if protein[i] in pl[j][1]:
#             if i == 0 and '_C' in pl[j][1]:
#                 ATOM.extend([
#                     [pl[j][0], pl[j][1], 0, 0, 0, pl[j][2], float(pl[j][3]),
#                      pl[j][4],
#                      float(pl[j][5]),
#                      0 + rad() * float(pl[j][5]) ** 0.333, 0 + rad() * float(pl[j][5]) ** 0.333,
#                      0 + rad() * float(pl[j][5]) ** 0.333, k + 1]])
#                 k+=1
#             if '_C' in pl[j][1] and i == 1:
#                 z = (-1) ** 0.5
#                 x = y = 2
#                 while isinstance(z, complex):  # 判断z是否为复数
#                     x = random.uniform(0, 3.8)
#                     y = random.uniform(-1.8, 1.8)
#                     z = (3.8 ** 2 - x ** 2 - y ** 2) ** 0.5
#                 # print(x, y, z)
#
#                 ATOM.extend([
#                     [pl[j][0], pl[j][1], ATOM[0][2] + x, ATOM[0][3] + y, ATOM[0][4] + z, pl[j][2],
#                      float(pl[j][3]),
#                      pl[j][4],
#                      float(pl[j][5]),
#                      ATOM[0][2] + x + rad() * float(pl[j][5]) ** 0.333,
#                      ATOM[0][3] + y + rad() * float(pl[j][5]) ** 0.333,
#                      ATOM[0][4] + z + rad() * float(pl[j][5]) ** 0.333, k + 1]])
#             if '_C' in pl[j][1] and i > 1:
#                 z = (-1) ** 0.5
#                 x = y = 2
#                 while isinstance(z, complex):  # 判断z是否为复数
#                     x = random.uniform(0, 3.8)
#                     y = random.uniform(-1.8, 1.8)
#                     z = (3.8 ** 2 - x ** 2 - y ** 2) ** 0.5
#                 # print(x, y, z)
#
#                 ATOM.extend([
#                     [pl[j][0], pl[j][1], ATOM[i - 2][2] + x, ATOM[i - 2][3] + y, ATOM[i - 2][4] + z, pl[j][2],
#                      float(pl[j][3]),
#                      pl[j][4],
#                      float(pl[j][5]),
#                      ATOM[i - 2][2] + x + rad() * float(pl[j][5]) ** 0.333,
#                      ATOM[i - 2][3] + y + rad() * float(pl[j][5]) ** 0.333,
#                      ATOM[i - 2][4] + z + rad() * float(pl[j][5]) ** 0.333, k + 1]])
#                 k += 1
#             if '_R' in pl[j][1]:
#                 # print(k)
#                 ATOM.extend([
#                     [pl[j][0], pl[j][1], ATOM[i - 1][2] + float(pl[j][3]), ATOM[i - 1][3] + float(pl[j][3]),
#                      ATOM[i - 1][4] + float(pl[j][3]), pl[j][2], float(pl[j][3]),
#                      pl[j][4],
#                      float(pl[j][5]),
#                      ATOM[i - 1][2] + float(pl[j][3]) + rad() * float(pl[j][5]) ** 0.333,
#                      ATOM[i - 1][3] + float(pl[j][3]) + rad() * float(pl[j][5]) ** 0.333,
#                      ATOM[i - 1][4] + float(pl[j][3]) + rad() * float(pl[j][5]) ** 0.333, k + 1]])
#                 k += 1
# print(ATOM)
# 以下为错误生成
# for i in range(len(protein)):
#     for j in range(len(pl)):
#         if protein[i] in pl[j][1]:
#             if '_C' in pl[j][1]:
#                 dis_two = 0
#                 if i == 0:
#                     cx = cy = cz = 0
#                 else:
#                     cx = random.uniform(0, 3.8)
#                     cy = random.uniform(-3.8, 3.8)
#                     cz = (3.8 ** 2 - cx ** 2 - cy ** 2) ** 0.5
#
#
#
#
#                 if i == 0:# 第一个原子在原点
#                     ATOM.extend([
#                         [pl[j][0], pl[j][1],
#                          0, 0, 0,  # 质心坐标
#                          pl[j][2], pl[j][3], pl[j][4], pl[j][5],
#                          rad() * (float(pl[j][5]) ** 0.333), rad() * (float(pl[j][5]) ** 0.333),
#                          rad() * (float(pl[j][5]) ** 0.333),  # 电荷坐标
#                          k + 1]])
#                 else :
#                     ATOM.extend([
#                         [pl[j][0], pl[j][1],
#                          ATOM[i - 2][2] + cx, ATOM[i - 2][3] + cy, ATOM[i - 2][4] + cz,  # 质心坐标
#                          pl[j][2], pl[j][3], pl[j][4], pl[j][5],
#                          cx + rad() * (float(pl[j][5]) ** 0.333), cy + rad() * (float(pl[j][5]) ** 0.333),
#                          cz + rad() * (float(pl[j][5]) ** 0.333),  # 电荷坐标
#                          k + 1]])
#
#             if '_R' in pl[j][1]:
#                 radry = rad()
#                 radrz = rad()
#                 ATOM.extend([
#                     [pl[j][0], pl[j][1],
#                      ATOM[i - 1][2], ATOM[i - 1][3] + radry * (float(pl[j][3])) ** 0.5, # 质心坐标
#                      ATOM[i - 1][4] + radrz * (float(pl[j][3])) ** 0.5,
#                      pl[j][2], pl[j][3], pl[j][4], pl[j][5],
#                      ATOM[i - 1][9] + float(pl[j][3]), ATOM[i - 1][10] + radry() * float(pl[j][3]),# 电荷坐标
#                      ATOM[i - 1][11] + radrz * float(pl[j][3]),
#                      k + 1]])
#             k += 1
# loadtest = open("loadtest.txt", 'w')
# for j in range(ATOM):
#     loadtest.write(
#         'ATOM' + '    ' + "{:>3d}".format(int(ATOM[j][12])) + '  ' + "{:<3s}".format(
#             str(ATOM[j][1])) + ' ' +
#         ATOM[j][0] + ' ' + 'A' + "{:>3d}".format(
#             math.ceil(int(ATOM[j][12]) / 2)) + '     ' + "{:>8.3f}".format(
#             float(ATOM[j][2])) + "{:>8.3f}".format(float(ATOM[j][3])) + "{:>8.3f}".format(
#             float(ATOM[j][4])) + '  ' + '1.00' + "{:>6.2f}".format(
#             134.34) + '           ' + ATOM[j][1][0] + ATOM[j][1][2] + '\n')

# 生成螺旋线
# for i in range(len(protein)):
#     j = 0
#     for j in range(len(pl)):
#         if protein[i] in pl[j][1]:
#             if '_C' in pl[j][1]:
#                 ATOM.extend([
#                     [pl[j][0], pl[j][1], i * 2, 16 * math.cos(i), 16 * math.sin(i), pl[j][2], pl[j][3], pl[j][4],
#                      pl[j][5],
#                      i * 2 + float(pl[j][5]), 16 * math.cos(i), 16 * math.sin(i), k + 1]])
#             if '_R' in pl[j][1]:
#                 ATOM.extend([
#                     [pl[j][0], pl[j][1], i * 2, 16 * math.cos(i) + 3.9, 16 * math.sin(i), pl[j][2], pl[j][3], pl[j][4],
#                      pl[j][5],
#                      i * 2 + float(pl[j][5]), 16 * math.cos(i) + 3.9, 16 * math.sin(i), k + 1]])
#             k += 1


# 在xoz平面上，延x轴震荡向正方向排序。 由于本程序粗粒力无垂直力，会扁平化扩散
# for i in range(len(protein)):
#     j = 0
#     for j in range(len(pl)):
#         if protein[i] in pl[j][1]:
#             ATOM.extend([
#                 [pl[j][0], pl[j][1], i * 2, 16 * math.cos(i), 16 * math.sin(i), pl[j][2], pl[j][3], pl[j][4],
#                  pl[j][5], i * 2 + float(pl[j][5]), 16 * math.cos(i), 16 * math.sin(i), k + 1]])
#             k += 1
#         j += 1
#     i += 1
#
# for i in range(len(ATOM)):
#     print(str(ATOM[i]))

# for i in range(len(ATOM)):
#     culi.write(str(ATOM[i]))
# print(ATOM)
for i in range(len(ATOM)):
    for j in range(len(ATOM[i])):
        culi.write(str(ATOM[i][j]))
        culi.write(' ')
    culi.write('\n')
culi.close()

connect = []
for i in range(len(ATOM)):
    if i == 0:
        connect.extend([[1, 2, 3]])
    elif i == len(ATOM) - 2:
        connect.extend([[i + 1, i - 1, i + 2]])
    elif (i + 1) % 2 == 0:
        connect.extend([[i + 1, i]])
    else:
        connect.extend([[i + 1, i - 1, i + 2, i + 3]])
# CON.write(str(connect))
for i in range(len(connect)):
    CON.write('CONECT' + '  ')
    for j in range(len(connect[i])):
        if connect[i][j] != '':
            CON.write("{:>3d}".format(int(connect[i][j])))
            CON.write('  ')
        else:
            break
    CON.write('\n')
CON.close()

for j in range(len(ATOM)):
    pdb.write(
        'ATOM' + '    ' + "{:>3d}".format(int(ATOM[j][12])) + '  ' + "{:<3s}".format(str(ATOM[j][1])) + ' ' +
        ATOM[j][0] + ' ' + 'A' + "{:>3d}".format(
            math.ceil(int(ATOM[j][12]) / 2)) + '     ' + "{:>8.3f}".format(
            float(ATOM[j][2])) + "{:>8.3f}".format(float(ATOM[j][3])) + "{:>8.3f}".format(
            float(ATOM[j][4])) + '  ' + '1.00' + "{:>6.2f}".format(
            134.34) + '           ' + ATOM[j][1][0] + ATOM[j][1][2] + '\n')
pdb.write('TER' + '\n')
for i in range(len(connect)):  # 向pdb文件写入connect
    pdb.write('CONECT' + '  ')
    for j in range(len(connect[i])):
        if connect[i][j] != '':
            pdb.write("{:>3d}".format(int(connect[i][j])))
            pdb.write('  ')
        else:
            break
    pdb.write('\n')
pdb.write('END')
pdb.close()

zero = np.array([0.0, 0.0, 0.0]).astype(np.float64)  # 定义零向量
low_energy = 0
# ------------------  低能量陷阱 -----------------------#
for cir in range(1):
    # if low_energy != 0:
    #     ATOM[0][2] = 5 * (-1) ** low_energy + ATOM[0][2]
    #

    # ------------------  低能量陷阱 -----------------------#
    angle_U = 0
    bound_U = 0
    vdW_U = 0
    dih_U = 0
    ele_U = 0
    # angle = zero
    # bound = zero

    m = 1  # 行进倍率
    step = 0

    BTOM = deepcopy(ATOM)
    start = time.process_time()
    while m:
        energy = 0
        U_x = 0
        ATOM = deepcopy(BTOM)

        dis_list = []  # 计算全局distance
        dis_list_x = []
        for j in range(len(ATOM)):
            for k in range(len(ATOM)):
                dis_list_x.extend([distance(ATOM[j], ATOM[k])])
            dis_list.extend([dis_list_x])
        # for j in range(len(dis_list)):
        #     print(str(dis_list[j]))
        #     print('\n\n')

        # for j in range(len(dis_list)):
        #     print(dis_list[j])
        dis_e_list = []  # 计算全局distance_e
        dis_e_list_x = []
        for j in range(len(ATOM)):
            for k in range(len(ATOM)):
                dis_e_list_x.extend([distance(ATOM[j], ATOM[k])])
            dis_e_list.extend([dis_e_list_x])

        for i in range(len(ATOM)):

            angle_U = 0
            bound_U = 0
            vdW_U = 0
            dih_U = 0
            ele_U = 0

            bound = deepcopy(zero)
            angle = deepcopy(zero)
            vdw = deepcopy(zero)
            dih = deepcopy(zero)
            ele = deepcopy(zero)

            force_all = deepcopy(zero)

            # 键能项
            for j in range(len(connect[i])):
                if j != 0 and connect[i][j] != '':
                    # print(connect[i][j] - 1)
                    # print(j)
                    number = connect[i][j] - 1  # 第二个原子序号
                    # print(i+1,number+1)
                    # print(i,j)
                    bound_U += U_bound(ATOM[i], ATOM[number])
                    # print(bound_U)
                    bound = force_bound(ATOM[i], ATOM[number]) + bound  # 注：向量计算不可用形如 a += 1 累加

            # 键角项
            if '_C' in ATOM[i][1]:
                if i == 0:
                    angle_U = U_angle_fix(ATOM[6], ATOM[7], ATOM[i], ATOM[i + 2], ATOM[i + 3])
                    angle = force_angle_fix(ATOM[6], ATOM[7], ATOM[i], ATOM[i + 2], ATOM[i + 3])
                if i == len(ATOM) - 2:
                    angle_U = U_angle_fix(ATOM[i - 2], ATOM[i - 1], ATOM[i], ATOM[6], ATOM[7])
                    angle = force_angle_fix(ATOM[i - 2], ATOM[i - 1], ATOM[i], ATOM[6], ATOM[7])
                else:
                    angle_U = U_angle_fix(ATOM[i - 2], ATOM[i - 1], ATOM[i], ATOM[i + 2], ATOM[i + 3])
                    angle = force_angle_fix(ATOM[i - 2], ATOM[i - 1], ATOM[i], ATOM[i + 2], ATOM[i + 3])

            # 以下为作用于侧链原子的键角项
            # if '_R' in ATOM[i][1] and ATOM[i][1][0] != 'G':  # 无关于Gly的高斯参数，已集成于函数
            #     if i == 1:  # 由于函数中已经对极端情况进行判断，第一个实参随便赋值
            #         angle_U = U_angle(ATOM[6], ATOM[i - 1], ATOM[i], ATOM[i + 1])
            #         angle = force_angle(ATOM[6], ATOM[i - 1], ATOM[i], ATOM[i + 1])
            #     if i == len(ATOM) - 1:  # 末端极端情况，第四个实参随便赋值
            #         angle_U = U_angle(ATOM[i - 3], ATOM[i - 1], ATOM[i], ATOM[6])
            #         angle = force_angle(ATOM[i - 3], ATOM[i - 1], ATOM[i], ATOM[6])
            #     else:
            #         angle_U = U_angle(ATOM[i - 3], ATOM[i - 1], ATOM[i], ATOM[i + 1])
            #         angle = force_angle(ATOM[i - 3], ATOM[i - 1], ATOM[i], ATOM[i + 1])

            # 注释部分为通过connect遍历法，与注释的键角函数配合
            # a = connect[i][0]  # a为被作用的atom编号
            # for j in range(len(connect[i])):
            #     if connect[i][j] != connect[i][0]:
            #         b = connect[i][j]  # b为做夹角交点的atom编号
            #         for k in range(len(connect[b - 1])):
            #             if connect[b - 1][k] != a and connect[b - 1][k] != b:
            #                 c = connect[b - 1][k]  # c为另一个atom编号
            #                 # print(a,b,c)
            #                 angle_U += U_angle(ATOM[b - 1], ATOM[a - 1], ATOM[c - 1])
            #                 angle = force_angle(ATOM[b - 1], ATOM[a - 1],
            #                                     ATOM[c - 1]) + angle  # 虚拟键角力  交点为atom1，作用于atom2，方向范围 作用强度为9

            # 范德华势能项     对大于20A的不考虑范德华势能
            for j in range(len(ATOM)):
                if dis_list[i][j] <= 20 and i != j:
                    vdW_U += U_vdW(ATOM[i], ATOM[j])
                    vdw = force_vdW(ATOM[i], ATOM[j]) + vdw

            # 二面角项 只考虑四个粒子中最后一个受力

            if i > 2 and '_C' in ATOM[i][1]:  # 只考虑a C的二面角能
                # 四个a C
                if i > 5:
                    dih_U += U_dih(ATOM[i - 6], ATOM[i - 4], ATOM[i - 2], ATOM[i])
                    dih = force_dih(ATOM[i - 6], ATOM[i - 4], ATOM[i - 2], ATOM[i]) + dih
                # 一个R
                dih_U += U_dih(ATOM[i - 3], ATOM[i - 4], ATOM[i - 2], ATOM[i])
                dih = force_dih(ATOM[i - 3], ATOM[i - 4], ATOM[i - 2], ATOM[i]) + dih

                # if '_R' in ATOM[i][1]:
                #     # RCCR
                #     dih_U += U_dih(ATOM[i - 2], ATOM[i - 3], ATOM[i - 1], ATOM[i])
                #     dih = force_dih(ATOM[i - 6], ATOM[i - 4], ATOM[i - 2], ATOM[i]) + dih
                #     # RCCC
                #     if i > 3:
                #         dih_U += U_dih(ATOM[i - 5], ATOM[i - 3], ATOM[i - 1], ATOM[i])
                #         dih = force_dih(ATOM[i - 5], ATOM[i - 4], ATOM[i - 2], ATOM[i]) + dih

            # 电荷能项
            for j in range(len(ATOM)):
                if dis_e_list[i][j] < 20 and j != i:
                    ele_U += U_ele(ATOM[i], ATOM[j])
                    ele = force_ele(ATOM[i], ATOM[j]) + ele
                    # print(i,j,ele,dis_e_list[i][j],ele_U,float(np.linalg.norm(ele))) # ele 调试

            force_all = bound + angle + vdw + dih + ele
            U_all = angle_U + bound_U + vdW_U + dih_U + ele_U
            force_all = force_all.astype(np.float64)
            f = float(np.linalg.norm(force_all))
            if float(ATOM[i][5]) == 0:
                acc = 0
            else:
                acc = f / float(ATOM[i][5])
            x = 0.5 * acc * 4 * m  # 向该方向前进1倍  步距 2 fs
            if x > 0.4:
                x = 0.4

            U_x += x

            # --------------检测各个力的数量级-------------#
            # print('----------', i, '----------')
            # print('bound=', bound, float(np.linalg.norm(bound)))
            # print('angle=', angle, float(np.linalg.norm(angle)))
            # print('vdw=', vdw, float(np.linalg.norm(vdw)))
            # print('dih', dih, float(np.linalg.norm(dih)))
            # print('ele=', ele, float(np.linalg.norm(ele)))
            # print('m=', m)
            # print('x=', x)
            # print('\n')
            # --------------检测各个力的数量级-------------#

            # 更新质心坐标
            place_1 = np.array([float(ATOM[i][2]), float(ATOM[i][3]), float(ATOM[i][4])])  # 原位置
            place_2 = place_1 + ((force_all / f) * x)
            if i == 0:
                place_2 = place_1
            BTOM[i][2] = place_2[0]
            BTOM[i][3] = place_2[1]
            BTOM[i][4] = place_2[2]

            # 更新电荷坐标
            place_ele = np.array([float(ATOM[i][9]), float(ATOM[i][10]), float(ATOM[i][11])])
            if float(np.linalg.norm(ele)) < 10 ** -9:
                place_ele_2 = place_ele
            else:
                place_ele_2 = place_ele + float(ATOM[1][8]) * ele / float(np.linalg.norm(ele))
            BTOM[i][9] = float(place_ele_2[0])
            BTOM[i][10] = float(place_ele_2[1])
            BTOM[i][11] = float(place_ele_2[2])

            energy += U_all  # 位置错误重写

        average = U_x / len(ATOM)
        if average < 0.2: m = 1
        if average < 0.1: m = 0.5
        if average < 0.01: m = 0
        step += 1
        if m == 0 or step % 100 == 0:
            # load = np.loadtxt('F:/毕业设计库/待检rmsd_norm.txt', dtype=str)  # 计算蛋白
            dir1 = np.array([[float(rmsd_norm[len(rmsd_norm) - 1][6]) - float(rmsd_norm[0][6])],
                             [float(rmsd_norm[len(rmsd_norm) - 1][7]) - float(rmsd_norm[0][7])],
                             [float(rmsd_norm[len(rmsd_norm) - 1][8]) - float(rmsd_norm[0][8])]])
            dir2 = np.array([[float(ATOM[len(rmsd_norm) - 1][2]) - float(ATOM[0][2])],
                             [float(ATOM[len(rmsd_norm) - 1][3]) - float(ATOM[0][3])],
                             [float(ATOM[len(rmsd_norm) - 1][4]) - float(ATOM[0][4])]])
            # print(dir1)
            # print(dir2)
            dir1_unit = dir1  # / np.linalg.norm(dir1)
            dir2_unit = dir1  # / np.linalg.norm(dir2)
            a = [float(dir1_unit[0]) / float(dir2_unit[0]), float(dir1_unit[1]) / float(dir2_unit[1]),
                 float(dir1_unit[2]) / float(dir2_unit[2])]
            p = np.diag(a)
            # print(str(p))
            add = 0
            for k in range(len(rmsd_norm)):
                space1 = np.array([[float(rmsd_norm[k][6]) - float(rmsd_norm[0][6])],
                                   [float(rmsd_norm[k][7]) - float(rmsd_norm[0][7])],
                                   [float(rmsd_norm[k][8]) - float(rmsd_norm[0][8])]])
                space2 = np.array([[float(ATOM[k][2])], [float(ATOM[k][3])], [float(ATOM[k][4])]])
                # print(p)
                # print(space2)
                transfer = np.matmul(p, space2)
                # print(transfer)
                # print(space1)
                dis = space1 - transfer
                # print(dis)
                add += float(dis[0]) ** 2 + float(dis[1]) ** 2 + float(dis[2]) ** 2
                # break
            rmsd = (add / len(rmsd_norm)) ** 0.5
            print('step=', step, 'm=', m, 'Move average=', average, 'energy=', energy, 'RMSD=', rmsd)

            end = time.process_time()
            sec = end - start
            note.write(str(sec) + ' ' + str(step) + ' ' + str(m) + ' ' + str(average) + ' ' + str(energy) + ' ' + str(
                rmsd) + '\n')

            if m == 0:
                f = open("{}end_{}_{}.txt".format(path, step, sec), 'w')
                result = open("{}end_{}_{}.pdb".format(path, step, sec), 'w')
            else:
                f = open("{}{}_{}.txt".format(path, step, sec), 'w')
                result = open("{}{}_{}.pdb".format(path, step, sec), 'w')

            for j in range(len(BTOM)):
                result.write(
                    'ATOM' + '    ' + "{:>3d}".format(int(BTOM[j][12])) + '  ' + "{:<3s}".format(
                        str(BTOM[j][1])) + ' ' +
                    BTOM[j][0] + ' ' + 'A' + "{:>3d}".format(
                        math.ceil(int(BTOM[j][12]) / 2)) + '     ' + "{:>8.3f}".format(
                        float(BTOM[j][2])) + "{:>8.3f}".format(float(BTOM[j][3])) + "{:>8.3f}".format(
                        float(BTOM[j][4])) + '  ' + '1.00' + "{:>6.2f}".format(
                        134.34) + '           ' + BTOM[j][1][0] + BTOM[j][1][2] + '\n')
                f.write(str(BTOM[j]))
                # f.write(
                #     'ATOM' + '    ' + "{:>3d}".format(int(BTOM[j][12])) + '  ' + "{:<3s}".format(
                #         str(BTOM[j][1])) + ' ' +
                #     BTOM[j][0] + ' ' + 'A' + "{:>3d}".format(
                #         math.ceil(int(BTOM[j][12]) / 2)) + '     ' + "{:>8.3f}".format(
                #         float(BTOM[j][2])) + "{:>8.3f}".format(float(BTOM[j][3])) + "{:>8.3f}".format(
                #         float(BTOM[j][4])) + '  ' + '1.00' + "{:>6.2f}".format(
                #         134.34) + '           ' + BTOM[j][1][0] + BTOM[j][1][2] + '\n')

            result.write('TER')
            result.write('\n')
            for k in range(len(connect)):  # 向pdb文件写入connect
                result.write('CONECT' + '  ')
                for j in range(len(connect[k])):
                    if connect[k][j] != '':
                        result.write("{:>3d}".format(int(connect[k][j])))
                        result.write('  ')
                    else:
                        break
                result.write('\n')
            result.write('END')
            f.close()
            result.close()
        # -------------------当前能量测试---------------------#
        # 加上后循环一次计算
        # print('step=', step, 'average=', average, 'energy=', energy)
        # if step == 50:
        #     break
    # -------------------当前能量测试---------------------#

note.close()
# #成键
#     if i>1:
#         k = i - 2
#         chengjian = []
#         while(k <= i + 2):
#             if ATOM[i][0] in connect[k]:
#                 chengjian.extend(connect[k][0])
#         k = k + 1
#     elif i == 0:
#         chengjian = [2,3]
#     else：
#         chengjian = [0]
#

# if i+1 in connect


# print(pl)
# b = pl[['摩尔质量']].values
# print(b)
# print(ATOM[2][1])
