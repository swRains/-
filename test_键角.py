import math

import numpy as np

# from fx import *
#
# w = open("F:/毕业设计库/angle_jj.pdb",'w')
#
# Atom = np.loadtxt('F:/毕业设计库/4dcz_ATOM_10.pdb',dtype= str)
# angle = angle_emj(Atom[1], Atom[2], Atom[3],Atom[4])
# print(angle)

BTOM = np.loadtxt('F:/毕业设计库/粗粒化后.txt', dtype=str)
result = open("F:/毕业设计库/write.pdb", 'w')
for j in range(len(BTOM)):
    result.write(
        'ATOM' + '    ' + "{:>3d}".format(int(BTOM[j][12])) + '  ' + "{:<3s}".format(str(BTOM[j][1])) +  ' ' + BTOM[j][0] + ' ' + 'A' + "{:>3d}".format(
            math.ceil(int(BTOM[j][12]) / 2)) + '     ' + "{:>8.3f}".format(
            float(BTOM[j][2])) + "{:>8.3f}".format(float(BTOM[j][3])) + "{:>8.3f}".format(
            float(BTOM[j][4])) + '  ' + '1.00' + "{:>6.2f}".format(
            134.34) + '           ' + BTOM[j][1][0] + BTOM[j][1][2] + '\n')
