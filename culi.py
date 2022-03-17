# 本代码为生成粗粒的程序
import numpy as np

w = open("F:/毕业设计库/culi.txt", 'w')
pl = np.loadtxt("F:/毕业设计库/protein_ele.txt",dtype = str)
Atom = np.loadtxt('F:/毕业设计库/culi_B.pdb', dtype=str)
j = 1
k = 1
for i in range(len(Atom)):
    if Atom[i][2] == 'CA':
        for j in range(len(pl)):
            if pl[j][0] == Atom[i][3] and '_C' in pl[j][1]:
                w.write(
                    Atom[i][0] + '    ' + "{:>3d}".format(k) + '  ' + "{:<3s}".format(pl[j][1]) + ' ' + Atom[i][
                        3] + ' ' + Atom[i][4] + ' ' + Atom[i][5] + '      ' + "{:>6.3f}".format(
                        float(Atom[i][6])) + '  ' + "{:>6.3f}".format(
                        float(Atom[i][7])) + '  ' + "{:>6.3f}".format(
                        float(Atom[i][8])) + '  ' + '1.00' + "{:>6.2f}".format(
                        float(134.34)) + '           ' + pl[j][1][0]+pl[j][1][2] + '\n')
        k += 1
    if Atom[i][3] == 'GLY':
        for j in range(len(pl)):
            if pl[j][0] == Atom[i][3] and '_C' in pl[j][1]:
                w.write(
                    Atom[i][0] + '    ' + "{:>3d}".format(k) + '  ' + "{:<3s}".format('G_R') + ' ' + Atom[i][
                        3] + ' ' + Atom[i][4] + ' ' + Atom[i][5] + '      ' + "{:>6.3f}".format(
                        float(Atom[i][6])) + '  ' + "{:>6.3f}".format(
                        float(Atom[i][7])) + '  ' + "{:>6.3f}".format(
                        float(Atom[i][8])) + '  ' + '1.00' + "{:>6.2f}".format(
                        float(134.34)) + '           ' + pl[j][1][0]+'R' + '\n')
        k+=1
    if Atom[i][2] == 'CB':
        for j in range(len(pl)):
            if pl[j][0] == Atom[i][3] and '_R' in pl[j][1]:
                w.write(
                    Atom[i][0] + '    ' + "{:>3d}".format(k) + '  ' + "{:<3s}".format(pl[j][1]) + ' ' + Atom[i][
                        3] + ' ' + Atom[i][4] + ' ' + Atom[i][5] + '      ' + "{:>6.3f}".format(
                        float(Atom[i][6])) + '  ' + "{:>6.3f}".format(
                        float(Atom[i][7])) + '  ' + "{:>6.3f}".format(
                        float(Atom[i][8])) + '  ' + '1.00' + "{:>6.2f}".format(
                        float(134.34)) + '           ' + pl[j][1][0] + pl[j][1][2] + '\n')
        k += 1

print('end')
