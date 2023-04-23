import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    with open("dots.txt") as f:
        m = int(f.readline())
        n = int(f.readline())
        u = [list(map(float, f.readline().split())) for i in range(n)]

    h = 1 / (m - 1)
    tau = 1 / (n - 1)

    x = np.arange(0, 1 + h, h)
    t = np.arange(0, 1 + tau, tau)

    x, y = np.meshgrid(x, t)

    plt.contourf(x, y, u)
    plt.title("Contour lines")
    plt.colorbar()
    plt.show()
