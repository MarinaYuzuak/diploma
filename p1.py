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

plt.scatter(50 , 50, c = 'deeppink')
plt.scatter(50 , 350, c = 'deeppink')

plt.scatter(300 , 200, c = 'black')
plt.scatter(850 , 200, c = 'black')

plt.grid()
plt.show()