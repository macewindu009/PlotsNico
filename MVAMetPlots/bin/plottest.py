import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 2 * np.pi, 100)
y1 = np.sin(x)
y2 = np.cos(3 * x)
plt.figure(0)
plt.plot(x,y1)
plt.figure(6)
plt.plot(x,y2)
plt.plot(x,y1)
plt.figure(0)
plt.show()
plt.savefig("test.png")
plt.figure(6)
plt.savefig("test1.png")
