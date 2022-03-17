from fx import *

file = open('F:/毕业设计库/culi_B.pdb', 'r')
result = open('F:/毕业设计库/angel.pdb', 'w')
list = []
line = file.readlines()
for i in range(len(line)):
    list.extend([line[i].split()])
# print(list[0])
# for i in range(len(list)):
#     print(list[i])
angle = []
for i in range(len(list)):
    if i == 0:
        angle.extend(
            [[list[i + 1][3], list[i][3], list[i + 2][3], 'R-C-C', angle_jj(list[i], list[i + 1], list[i + 2])]])
    elif i == len(list) - 2:
        angle.extend(
            [[list[i - 2][3], list[i][3], list[i + 1][3], 'C-C-R', angle_jj(list[i], list[i - 2], list[i + 1])]])
    elif list[i][2] == 'CA':
        # print(i)
        angle.extend(
            [[list[i - 2][3], list[i][3], list[i + 1][3], 'C-C-R', angle_jj(list[i], list[i - 2], list[i + 1])]])
        angle.extend(
            [[list[i + 1][3], list[i][3], list[i + 2][3], 'R-C-C', angle_jj(list[i], list[i + 1], list[i + 2])]])
        angle.extend(
            [[list[i - 2][3], list[i][3], list[i + 2][3], 'C-C-C', angle_jj(list[i], list[i - 2], list[i + 2])]])

    # elif list[i][3] == 'GLY' and list[i][2] == 'CA':
    #     angle.extend([list[i - 2][3], list[i][3], list[i + 1][3], angle_jj(list[i], list[i - 2], list[i + 1])])

for i in range(len(angle)):
    print(angle[i])
    result.write(angle[i][0] + '  ' + angle[i][1] + '  ' + angle[i][2] + '  ' + angle[i][3] + '  ' + "{:<10.6f}".format(
        angle[i][4]) + '  ' + '\n')

print(angle[7][4])

# result.write(angle[i][0] + '  ' + angle[i][1] + angle[i][2] + angle[i][3] + "{:<10.6f}".format(angle[i][4]))
