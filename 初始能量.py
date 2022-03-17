from fx_fix import *
ATOM = []
k = 1
ortxt = open('F:/毕业设计库/or.txt',"w")
origion = np.loadtxt('F:/毕业设计库/origion.pdb',dtype = str)
pl = np.loadtxt('F:/毕业设计库/protein_ele.txt', dtype=str)  # protein list
connect = []

for i in range(len(origion)):
    for j in range(len(pl)):
        if pl[j][1] == origion[i][2]:
            ATOM.extend([
                    [pl[j][0], pl[j][1], float(origion[i][6]), float(origion[i][7]), float(origion[i][8]), pl[j][2],
                     float(pl[j][3]),
                     pl[j][4],
                     float(pl[j][5]),
                     float(origion[i][6]), float(origion[i][7]), float(origion[i][8]), k]])
            k += 1
        if origion[i][2] == 'G_R':
            ATOM.extend([
                [pl[j][0], pl[j][1], float(origion[i][6]), float(origion[i][7]), float(origion[i][8]), pl[j][2],
                 float(pl[j][3]),
                 pl[j][4],
                 float(pl[j][5]),
                 float(origion[i][6]), float(origion[i][7]), float(origion[i][8]), k]])
            k += 1
ortxt.write(str(ATOM))
for i in range(len(ATOM)):
    if i == 0:
        connect.extend([[1, 2, 3]])
    elif i == len(ATOM) - 2:
        connect.extend([[i + 1, i - 1, i + 2]])
    elif (i + 1) % 2 == 0:
        connect.extend([[i + 1, i]])
    else:
        connect.extend([[i + 1, i - 1, i + 2, i + 3]])
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
energy = 0
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


    # 键角项
    if '_C' in ATOM[i][1]:
        if i == 0:
            angle_U = U_angle_fix(ATOM[6], ATOM[7], ATOM[i], ATOM[i + 2], ATOM[i + 3])

        if i == len(ATOM) - 2:
            angle_U = U_angle_fix(ATOM[i - 2], ATOM[i - 1], ATOM[i], ATOM[6], ATOM[7])

        else:
            angle_U = U_angle_fix(ATOM[i - 2], ATOM[i - 1], ATOM[i], ATOM[i + 2], ATOM[i + 3])


    # 范德华势能项     对大于20A的不考虑范德华势能
    for j in range(len(ATOM)):
        if dis_list[i][j] <= 20 and i != j:
            vdW_U += U_vdW(ATOM[i], ATOM[j])


    # 二面角项 只考虑四个粒子中最后一个受力

    if i > 2 and '_C' in ATOM[i][1]:  # 只考虑a C的二面角能
        # 四个a C
        if i > 5:
            dih_U += U_dih(ATOM[i - 6], ATOM[i - 4], ATOM[i - 2], ATOM[i])

        # 一个R
        dih_U += U_dih(ATOM[i - 3], ATOM[i - 4], ATOM[i - 2], ATOM[i])


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

            # print(i,j,ele,dis_e_list[i][j],ele_U,float(np.linalg.norm(ele))) # ele 调试


    U_all = angle_U + bound_U + vdW_U + dih_U + ele_U
    energy += U_all
print(energy)