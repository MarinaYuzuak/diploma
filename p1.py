import numpy as np
import matplotlib.pyplot as plt

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

plt.scatter(300 , 300, c = 'deeppink')
plt.scatter(300 , 700, c = 'deeppink')

plt.scatter(500 , 600, c = 'orange')
plt.scatter(600 , 600, c = 'orange')

plt.grid()
plt.show()