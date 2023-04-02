import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    x = np.arange(-1.9, 2.0, 0.01)
    # y = np.arange(-1.0, 2.0, 0.01)

    y1 = x ** 2
    y2 = np.log(x + 2)
    plt.plot(x, y1, 'r', label='x^2')
    plt.plot(x, y2, 'b', label='log(x+2)')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('y(x)')
    plt.grid()
    plt.legend()
    plt.show()

    x_s = np.arange(-2.0, 2.0, 0.01)
    y_s = np.arange(-1.0, 2.0, 0.01)

    x, y = np.meshgrid(x_s, y_s)

    z = (y - x ** 2) ** 2 + (x - np.exp(y) + 2) ** 2

    levels = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 5, 10, 20]

    plt.contourf(x, y, z, levels=levels)
    plt.title("Contour lines")
    plt.colorbar()
    plt.show()
