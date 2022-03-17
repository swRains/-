# import matplotlib.pyplot as plt
# import numpy as np
# import math
#
# x = np.arange(8,20,0.1)
# y = []
# z = []
# for i in x:
#     y_1 = 16 * math.cos(i)
#     y.append(y_1)
#     z_1 = 16 * math.sin(i)
#     z.append(z_1)
#
# plt.plot(x, y, label="vdW")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.ylabel("z")
# plt.ylim(0, 1)
# plt.legend()
# plt.show()

# import numpy as np
# import math
# import matplotlib.pyplot as plt
# x = np.arange(-10, 10, 0.1)
# y = []
# for t in x:
#     y_1 = 1 / (1 + math.exp(-t))
#     y.append(y_1)
# plt.plot(x, y, label="sigmoid")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.ylim(0, 1)
# plt.legend()
# plt.show()
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#定义图像和三维格式坐标轴
fig=plt.figure()
ax1 = Axes3D(fig)
import numpy as np
x = 3.8*np.linspace(0,100,100)
y = 0*x
z = 3 * np.sin(x/2)

xd = 3.8*np.linspace(0,100,100)
yd = 0*x
zd = 3 * np.sin(x/2)
ax1.scatter3D(xd,yd,zd, cmap='Blues')  #绘制散点图
ax1.plot3D(x,y,z,'gray')    #绘制空间曲线
plt.show()