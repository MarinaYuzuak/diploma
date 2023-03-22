import numpy as np
import matplotlib.pyplot as plt

data = []
with open("z.txt") as f:
    for line in f:
        data.append([float(x) for x in line.split()])

for d in data:
    plt.axhline(y=d, xmin=0, xmax=1)


plt.grid()
plt.show()
