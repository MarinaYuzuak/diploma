import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

fig, ax = plt.subplots()

x_data = []
with open("x.txt") as f:
    for line in f:
        x_data.append([float(x) for x in line.split()])

y_data = []
with open("y.txt") as f:
    for line in f:
        y_data.append([float(x) for x in line.split()])

for d in x_data:
    plt.axvline(x = d, ymin=0, ymax=1)


for d in y_data:
    plt.axhline(y=d, xmin=0, xmax=1)

ax.add_patch (Rectangle((1, 1), 2, 6))

plt.scatter(-500 , -1500, c = 'deeppink', s=30)
plt.scatter(1000 , -1500, c = 'deeppink', s=30)

#plt.scatter(-500 , -1580, c = 'orange', s=30)
#plt.scatter(1000 , -1580, c = 'orange', s=30)

#plt.scatter(-500 , -1710, c = 'orange', s=30)
#plt.scatter(1000 , -1710, c = 'orange', s=30)

plt.scatter(-500 , -1630, c = 'orange', s=30)
plt.scatter(1000 , -1630, c = 'orange', s=30)

plt.scatter(-500 , -1770, c = 'orange', s=30)
plt.scatter(1000 , -1770, c = 'orange', s=30)


plt.show()

