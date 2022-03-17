from fx_fix import *

pl = np.loadtxt('F:/毕业设计库/protein.txt', dtype=str)  # protein list
loadtest = open("loadtest.txt", 'w')
# protein = ['K_', 'Q_', 'E_', 'Q_', 'P_', 'E_', 'I_', 'N_', 'L_', 'D_', 'H_', 'V_', 'V_', 'E_', 'Q_', 'T_', 'I_', 'K_',
#            'K_', 'V_', 'Q_', 'Q_', 'N_', 'Q_', 'N_', 'Q_', 'N_', 'K_', 'D_', 'P_', 'D_', 'E_', 'L_', 'R_', 'S_', 'K_',
#            'V_', 'P_', 'G_', 'E_', 'V_', 'T_', 'A_', 'S_', 'D_', 'W_', 'E_', 'A_', 'L_', 'V_', 'G_', 'D_', 'T_', 'R_',
#            'Y_', 'G_', 'Y_', 'F_', 'D_', 'E_', 'T_', 'G_', 'D_', 'W_', 'S_', 'W_', 'K_', 'G_', 'Y_', 'F_', 'D_', 'E_',
#            'Q_', 'G_', 'K_', 'W_', 'V_', 'W_', 'N_', 'E_', 'P_', 'V_', 'D_', 'S_', 'L_', 'E_', 'H_', 'H_', 'H_', 'H_',
#            'H_', 'H_']
#
# # 以下为pdb中的蛋白序列 共计55个氨基酸 为26-80的切片
# protein = protein[25:79]  # 具体如下
protein = ['Q_', 'N_', 'K_', 'D_', 'P_', 'D_', 'E_', 'L_', 'R_', 'S_', 'K_', 'V_', 'P_', 'G_', 'E_', 'V_', 'T_', 'A_',
           'S_', 'D_', 'W_', 'E_', 'A_', 'L_', 'V_', 'G_', 'D_', 'T_', 'R_', 'Y_', 'G_', 'Y_', 'F_', 'D_', 'E_', 'T_',
           'G_', 'D_', 'W_', 'S_', 'W_', 'K_', 'G_', 'Y_', 'F_', 'D_', 'E_', 'Q_', 'G_', 'K_', 'W_', 'V_', 'W_', 'N_',
           'E_']
tzxl = []  # 拓展向量
dis_two = 0
# for i in range(len(protein)):
#     while dis_two > 3.6 and dis_two < 4.0:
#         if i == 0:
#             cx = cy = cz = 0
#         else:
#             cx = random.uniform(0, 3.8)
#             cy = random.uniform(-3.8, 3.8)
#             cz = random.uniform(-3.8, 3.8)
#             dis_two = (cx ** 2 + cy ** 2 + cz ** 2) ** 0.5
#         tzxl.append([[cx,cy,cz]])
# print(tzxl)

ATOM = []
# ATOM.extend([0])
k = 1
i = 0
x = y = z = 2
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
                    [pl[j][0], pl[j][1], ATOM[k-3][2] + x, ATOM[k-3][3] + y, ATOM[k-3][4] + z, pl[j][2],
                     float(pl[j][3]),
                     pl[j][4],
                     float(pl[j][5]),
                     ATOM[k-3][2] + x + rad() * float(pl[j][5]) ** 0.333,
                     ATOM[k-3][3] + y + rad() * float(pl[j][5]) ** 0.333,
                     ATOM[k-3][4] + z + rad() * float(pl[j][5]) ** 0.333, k]])

            ATOM.extend([
                [pl[j + 1][0], pl[j + 1][1], ATOM[k-1][2] + float(pl[j + 1][3]), ATOM[k-1][3] + float(pl[j + 1][3]),
                 ATOM[k-1][4] + float(pl[j + 1][3]), pl[j + 1][2], float(pl[j + 1][3]),
                 pl[j + 1][4],
                 float(pl[j + 1][5]),
                 ATOM[k-1][2] + float(pl[j + 1][3]) + rad() * float(pl[j + 1][5]) ** 0.333,
                 ATOM[k-1][3] + float(pl[j + 1][3]) + rad() * float(pl[j + 1][5]) ** 0.333,
                 ATOM[k-1][4] + float(pl[j + 1][3]) + rad() * float(pl[j + 1][5]) ** 0.333, k + 1]])
            k += 2

# for i in range(len(protein)):
#     j = 0
#     for j in range(len(pl)):
#         if protein[i] in pl[j][1]:
#             if i == 0:
#                 ATOM.extend([
#                     [pl[j][0], pl[j][1], 0, 0, 0, pl[j][2], float(pl[j][3]),
#                      pl[j][4],
#                      float(pl[j][5]),
#                      0 + rad() * float(pl[j][5]) ** 0.333, 0 + rad() * float(pl[j][5]) ** 0.333,
#                      0 + rad() * float(pl[j][5]) ** 0.333, k + 1]])
#             if '_C' in pl[j][1] and i >= 2:
#                 z = (-1) ** 0.5
#                 while isinstance(z, complex):# 判断z是否为复数
#                     x = random.uniform(0, 3.8)
#                     y = random.uniform(-1.8, 1.8)
#                     z = (3.8 ** 2 - x ** 2 - y ** 2) ** 0.5
#                 print(x, y, z)
#
#                 ATOM.extend([
#                     [pl[j][0], pl[j][1], ATOM[i - 2][2] + x, ATOM[i - 2][3] + y, ATOM[i - 2][4] + z, pl[j][2], float(pl[j][3]),
#                      pl[j][4],
#                      float(pl[j][5]),
#                      ATOM[i - 2][2] + x + rad() * float(pl[j][5]) ** 0.333, ATOM[i - 2][3] + y + rad() * float(pl[j][5]) ** 0.333,
#                      ATOM[i - 2][4] + z + rad() * float(pl[j][5]) ** 0.333, k + 1]])
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
#             k += 1
# for i in range(len(protein)):
#     for j in range(len(pl)):
#         # if i % 2 == 0:
#         #     cx = random.uniform(0, 3.8)
#         #     cy = random.uniform(-3.8, 3.8)
#         #     cz = (3.8 ** 2 - cx ** 2 - cy ** 2) ** 0.5
#         if i == 0:
#             if protein[i] in pl[j][1] and '_C' in pl[j]:  # 第一个原子在原点
#                 ATOM.extend([
#                     [pl[j][0], pl[j][1],
#                      0, 0, 0,  # 质心坐标
#                      pl[j][2], float(pl[j][3]), pl[j][4], float(pl[j][5]),
#                      rad() * (float(float(pl[j][5])) ** 0.333), rad() * (float(float(pl[j][5])) ** 0.333),
#                      rad() * (float(float(pl[j][5])) ** 0.333),  # 电荷坐标
#                      k + 1]])
#                 k += 1
#
#         elif protein[i] in pl[j][1] and '_C' in pl[j]:
#             cx = random.uniform(0, 3.8)
#             cy = random.uniform(-3.8, 3.8)
#             cz = (3.8 ** 2 - cx ** 2 - cy ** 2) ** 0.5
#             ATOM.extend([
#                 [pl[j][0], pl[j][1],
#                  ATOM[i - 2][2] + cx, ATOM[i - 2][3] + cy, ATOM[i - 2][4] + cz,  # 质心坐标
#                  pl[j][2], float(pl[j][3]), pl[j][4], float(pl[j][5]),
#                  cx + rad() * (float(float(pl[j][5])) ** 0.333), cy + rad() * (float(float(pl[j][5])) ** 0.333),
#                  cz + rad() * (float(float(pl[j][5])) ** 0.333),  # 电荷坐标
#                  k + 1]])
#             k += 1
#         if i % 2 != 0 and protein[i] in pl[j][1] and '_R' in pl[j]:
#             radr_y = rad()
#             radr_z = rad()
#             # print(ATOM)
#             ATOM.extend([
#                 [pl[j][0], pl[j][1],
#                  ATOM[i - 1][2], ATOM[i - 1][3] + radr_y * (float(float(pl[j][3]))) ** 0.5,  # 质心坐标
#                  ATOM[i - 1][4] + radr_z * (float(float(pl[j][3]))) ** 0.5,
#                  pl[j][2], float(pl[j][3]), pl[j][4], float(pl[j][5]),
#                  ATOM[i - 1][9] + float(float(pl[j][3])),
#                  ATOM[i - 1][10] + radr_z * float(float(pl[j][3])),  # 电荷坐标
#                  ATOM[i - 1][11] + radr_z * float(float(pl[j][3])),
#                  k + 1]])
#             k += 1

print(ATOM)
# print(1)
# for i in range(len(ATOM)):
#     print(ATOM[i])
#     loadtest.write(str(ATOM[i]) + '\n')
for j in range(len(ATOM)):
    # print(j,ATOM[j][4])
    loadtest.write(
        'ATOM' + '    ' + "{:>3d}".format(int(ATOM[j][12])) + '  ' + "{:<3s}".format(
            str(ATOM[j][1])) + ' ' +
        ATOM[j][0] + ' ' + 'A' + "{:>3d}".format(
            math.ceil(int(ATOM[j][12]) / 2)) + '     ' + "{:>8.3f}".format(
            float(ATOM[j][2])) + "{:>8.3f}".format(float(ATOM[j][3])) + "{:>8.3f}".format(
            float(ATOM[j][4]))
        + '  ' + '1.00' + "{:>6.2f}".format(
            134.34) + '           ' + ATOM[j][1][0] + ATOM[j][1][2] + '\n')

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
loadtest.write('TER' + '\n')
for i in range(len(connect)):  # 向pdb文件写入connect
    loadtest.write('CONNECT' + '  ')
    for j in range(len(connect[i])):
        if connect[i][j] != '':
            loadtest.write("{:>3d}".format(int(connect[i][j])))
            loadtest.write('  ')
        else:
            break
    loadtest.write('\n')
loadtest.write('END')
loadtest.close()
