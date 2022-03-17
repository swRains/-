import numpy as np

# import fx_fix
pdb = np.loadtxt('F:/毕业设计库/标准蛋白.txt', dtype=str)  # 标准蛋白
load = np.loadtxt('F:/毕业设计库/待检pdb.txt', dtype=str)  # 计算蛋白
dir1 = np.array([[float(pdb[len(pdb) - 1][6]) - float(pdb[0][6])], [float(pdb[len(pdb) - 1][7]) - float(pdb[0][7])],
                 [float(pdb[len(pdb) - 1][8]) - float(pdb[0][8])]])
dir2 = np.array([[float(load[len(pdb) - 1][6]) - float(load[0][6])], [float(load[len(pdb) - 1][7]) - float(load[0][7])],
                 [float(load[len(pdb) - 1][8]) - float(load[0][8])]])
# print(dir1)
# print(dir2)
dir1_unit = dir1 #/ np.linalg.norm(dir1)
dir2_unit = dir1 #/ np.linalg.norm(dir2)
a = [float(dir1_unit[0]) / float(dir2_unit[0]), float(dir1_unit[1]) / float(dir2_unit[1]),
     float(dir1_unit[2]) / float(dir2_unit[2])]
p = np.diag(a)
# print(str(p))
add = 0
for i in range(len(pdb)):
    space1 = np.array([[float(pdb[i][6]) - float(pdb[0][6])], [float(pdb[i][7]) - float(pdb[0][7])],
                       [float(pdb[i][8]) - float(pdb[0][8])]])
    space2 = np.array([[float(load[i][6])], [float(load[i][7])], [float(load[i][8])]])
    # print(p)
    # print(space2)
    transfer = np.matmul(p,space2)
    # print(transfer)
    # print(space1)
    dis = space1 - transfer
    # print(dis)
    add += float(dis[0]) ** 2 + float(dis[1]) ** 2 + float(dis[2]) ** 2
    # break
rmsd = (add / len(pdb)) ** 0.5
print('rmsd=' + str(rmsd))
