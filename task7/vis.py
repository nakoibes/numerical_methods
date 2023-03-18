import numpy as np
import matplotlib.pyplot as plt
x = np.arange(-2.0, 2.0, 0.01)
y = np.arange(-1.0, 1.0, 0.01)
y1 = x**2-1
x2 = np.tan(y)
plt.plot(x, y1, 'r', label='f(x, y) = x^2 - y - 1')
plt.plot(x2, y, 'b', label='g(x, y) = x - tg y')
plt.plot(1.396, 0.949, 'o')
plt.annotate('(1.396, 0.949)', (1.396 + 0.1, 0.949 + 0.1))

plt.xlabel('x')
plt.ylabel('y')
plt.title('y(x)')
plt.grid()
plt.legend()
plt.show()

x = np.arange(-2.0, 2.0, 0.01)
y = np.arange(-1.0, 1.0, 0.01)

xg, yg = np.meshgrid(x, y)

zg = (xg**2 - yg-1)**2 + (xg - np.tan(yg))**2

plt.contour(xg, yg, zg, levels = 200)
plt.title("Disperancy level's sum")
plt.colorbar()
plt.show()