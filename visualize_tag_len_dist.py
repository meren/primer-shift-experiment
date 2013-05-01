import sys
import matplotlib.pyplot as plt

f = open(sys.argv[1])

x = []
y = []

for line in f.readlines():
    _x, _y = line.strip().split()
    x.append(int(_x))
    y.append(int(_y))

plt.plot(x, y)
plt.xlim(xmin = min(x), xmax = max(x))
plt.ylim(ymax = max(y) * 1.1)
plt.xticks(range(min(x), max(x)))
plt.show()