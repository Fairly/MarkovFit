__author__ = 'Administrator'

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y, Z = axes3d.get_test_data(0.05)
# ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
# cset = ax.contour(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
# cset = ax.contour(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
# cset = ax.contour(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)
#
# ax.set_xlabel('X')
# ax.set_xlim(-40, 40)
# ax.set_ylabel('Y')
# ax.set_ylim(-40, 40)
# ax.set_zlabel('Z')
# ax.set_zlim(-100, 100)


def f(x, y):
    return np.sqrt(np.square(x)+np.square(y))


x = np.linspace(-1, 1, num=100)
y = np.linspace(-10, 10, num=100)
x, y = np.meshgrid(x, y)
z = f(x, y)

var = ax.plot_surface(x, y, z)
plt.show()