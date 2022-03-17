# encoding=utf-8
import numpy as np
import random
# import datetime
import time
import pandas as pd
#from fx_fix import *
import os
# f = open("F:/毕业设计库/4dcz_ATOM.pdb")
#w = open("F:/毕业设计库/cccc.txt", 'w')
# lines = f.readlines()
#Atom = np.loadtxt('F:/毕业设计库/cccc.txt', dtype=str)
# list = pd.read_csv('F:/毕业设计库/cccc.txt',index_col=0,sep=' ')
# #print(list)
# a =  list.loc['a1']
# b = a+1
# b = float(b)
# print(b)
#print(list.columns)
#print(list.iloc[0,'B-B-B-B'])
# data = pd.read_scv('F:/毕业设计库/4dcz_ATOM.pdb', sep=' ')
# print(data)
# i=0
# while(i<921):
#     j=i+1
#     while(j<=921):
#         dis = distance(Atom[i],Atom[j])#((float(Atom[i][6])-float(Atom[j][6]))**2+(float(Atom[i][7])-float(Atom[j][7]))**2+(float(Atom[i][8])-float(Atom[j][8]))**2)**0.5
#         if float(dis) <= 3:
#             #w.write(str(Atom[i][10])+' ---- '+str(Atom[j][10])+'      '+str(dis)+'\n')
#             print(dis)
#         #else:
#          #   print(i)
#         j+=1
#     print(i)
#     i+=1
# f.close()
# w.close()
# w.write("{:>3d}".format(Atom[2]))
# w.write(
#     Atom[0] + '     ' + "{:>2d}".format(int(Atom[1])) + '  ' + "{:<3s}".format(Atom[2]) + ' ' + Atom[3] + ' ' + Atom[
#         4] + ' ' + Atom[5] + '      ' + "{:>.3f}".format(float(Atom[6])) + '  ' + "{:>.3f}".format(
#         float(Atom[7])) + '  ' + "{:>.3f}".format(float(Atom[8])) + '  ' + '1.00' + "{:>.2f}".format(
#         float(Atom[9])) + '           ' + Atom[10])
#
# fx = force_vdW_dir(Atom[1], Atom[2])
# a = np.array([5,3,7])
# z = np.dot(a,fx)
# print(fx)
# print(z)
# x = angle_emj(atom1, atom2, atom3, atom4)
# i = 0
# while (i < 5):
#     U = dih_cccc[i] * math.exp(-((x - dih_cccc[i + 5]) / dih_cccc[i + 10]) ** 2)
#     i += 1
# test = [0]#错误
# b = int(test) + 1
# print(b)
#a = []
# a.extend([0])
# print(a[0])
# i = 0
# while i<5:
#     print(i)
#     i+= 1

# ccr = pd.read_csv('F:/毕业设计库/ccr.txt',index_col=0,sep=' ')
# print(ccr)
#
# rcc = pd.read_csv('F:/毕业设计库/rcc.txt',index_col=0,sep=' ')
# print(rcc)
#
# print(ccr.loc['a5', 'D'])


# a = 'C_R'
#
# print(a[0])

# a = [['adas','adwa','dawdawd'],
#      ['xmcnbwje','dajsgda','dajshgdj'],
#      ['sadw','uytq7wye','jahgkswi']]
# print(a[0][1])
# b = a[0][1][3]
# print(b)


# ccr = pd.read_csv('F:/毕业设计库/ccr.txt',index_col=0,sep=' ')
# print(ccr)
# record = 'CSJKHD'
# a = ccr[record[0]][0]
# print(a)
#
# for i in range(5):
#     f = open("{}.txt".format(i),'w')
#     f.write(str(i))

# start = time.perf_counter()
#
# end = time.perf_counter()
# time =  time.ctime()
# print(time,start,end,end-start)
# a = False
# i=0
# if i>10 and not a:
#     i +=1
#     print(i)
# if not a:
#     print(False)
# while 0:
#     print('no')
#     break
# while -0.001:
#     print('yes')
#     break
# print(time.localtime(time.time()))
# print(time.time())
# print(time.asctime( time.localtime(time.time()) ))
# print (time.strftime("%Y_%m_%d", time.localtime()) )
# localtime = time.strftime("%Y_%m_%d", time.localtime())
# path = 'F:/毕业设计库/测试数据/{:s}/'.format(localtime)
# os.makedirs(path)
# test = open('{}test.txt'.format(path), 'w')
# test.write('000')
# test.close()


# def mkdir(path):
#     folder = os.path.exists(path)
#
#     if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
#         os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径
#         print("---  新目录地址...  ---")
#         print("---  已创建新目录  ---")
#
#     else:
#         print("---  地址已存在，不再创建  ---")
#
# mkdir(path)
# for j in range(9):
#     a = [2, 3.7, 6, 9]
#     b = a.extend([5])
#     print(b)
#
#     asda = 3
# e = 1.6 * 10 ** -19
# U =6.02 * 10 ** 23 * e ** 2 * float(-0.063) * float(0.063) * 10 ** -3 / (4 * 3.1415926 * 8.85418 * 10 ** -12 * 10 ** -10 * 3)
# print(U)
# print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
# x1 = 1
# y1 = 2
# z1 = 3
# x2 = 4
# y2 = 5
# z2 = 6
# b = (x2 * z1 - z2 * x1) / (y2 * x1 - x2 * y1)
#
# a = -(b * y1 + z1) / x1
# c = 1
# print(a,b,c)

# f1 = np.array([1, 2, 3])
# f2 = f1 / np.linalg.norm(f1)
# print(f2[0],f2[1],f2[2])

# x1 = 1
# y1 = 2
# z1 = 4
#
# x2 = 2
# y2 = 4
# z2 = 7
#
# x3 = 3
# y3 = 6
# z3 = 11
#
# A = np.array([[x1, y1, z1],
#               [x2, y2, z2],
#               [x3, y3, z3]])
# B = np.array([[0], [0], [0]])
# if np.linalg.det(A) != 0:
#     x = B
#     print(x)
# elif np.linalg.det(A) == 0:
#     if np.linalg.matrix_rank(A) == 3:
#         x = np.linalg.solve(A, B)
#     if np.linalg.matrix_rank(A) == 2:
#     if np.linalg.matrix_rank(A) == 1:
#         x = 1
#     print(x)
# x1 = y1 = z1 = 0
# zero_forbid = np.array([0, 0, 0]).astype(np.float64)
# if x1 == 0: zero_forbid[0] = 1
# if y1 == 0: zero_forbid[1] = 1
# if z1 == 0: zero_forbid[2] = 1
#
# fx = zero_forbid / np.linalg.norm(zero_forbid)
# print(fx,type(fx))
#
# start = time.process_time()
# for i in range(100):
#     test = np.loadtxt('F:/毕业设计库/粗粒数据.txt',dtype = str)
# end = time.process_time()
# print (end - start)
# for i in range(50):
#     print(random.randint(-1,1))
# dis_two = 0.0
# z = (-1)**0.5
# while isinstance(z,complex):
#     x = y = z = 0
#     x = random.uniform(0, 3.8)
#     y = random.uniform(-3.8, 3.8)
#     z = (3.8 ** 2 - x ** 2 - y ** 2) ** 0.5
# print(x, y, z)

    # # if i == 0:
    # #     cx = cy = cz = 0
    # # else:
    #     cx = random.uniform(0, 3.8)
    #     cy = random.uniform(-3.8, 3.8)
    #     cz = random.uniform(-3.8, 3.8)
    #     dis_two = (cx ** 2 + cy ** 2 + cz ** 2) ** 0.5
    #     print(cx, cy, cz)

# print(cx,cy,cz)
# i=0
# while i<50:
#     i+=1
#     print(i)
# for i in range(100):
#     z = (-1) ** 0.5
#     x = y = 2
#     while isinstance(z, complex):  # 判断z是否为复数
#         x = random.uniform(0, 3.8)
#         y = random.uniform(-1.8, 1.8)
#         z = (3.8 ** 2 - x ** 2 - y ** 2) ** 0.5
#     print(x, y, z)

