import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

points = np.array([[400, 0, 0],
                      [400, 500, 0 ],
                      [0, 0, 0],
                      [0, 500, 0],
                      [400, 0, 600],
                      [400, 500, 600 ],
                      [0, 0, 600],
                      [0, 500, 600]])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
r = [400,600]
X, Y = np.meshgrid(r, r)
one = np.ones(4).reshape(2, 2)
ax.plot_wireframe(X,Y,600, alpha=0.5)
ax.plot_wireframe(X,Y,-600, alpha=0.5)
ax.plot_wireframe(X,-600,Y, alpha=0.5)
ax.plot_wireframe(X,600,Y, alpha=0.5)
ax.plot_wireframe(600,X,Y, alpha=0.5)
ax.plot_wireframe(-600,X,Y, alpha=0.5)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()