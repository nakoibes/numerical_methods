import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    with open("dots.txt") as f:
        k = int(f.readline())
        z1 = [[]] * k
        z2 = [[]] * k
        for i in range(k):
            z1[i] = list(map(int, f.readline().split()))
        for i in range(k):
            z2[i] = list(map(int, f.readline().split()))

    x_s = np.arange(-1.1, -0.1, 1.0 / k)
    y_s = np.arange(-0.15, 0.85, 1.0 / k)

    x, y = np.meshgrid(x_s, y_s)
    levels = [5, 10, 20, 50, 70, 200, 500, 2000]
    colors = ["red", "orange", "yellow", "green", "blue", "purple", "black"]
    plt.contourf(x, y, z1, levels=levels, colors=colors)
    plt.title("Root analysis")
    plt.colorbar()
    plt.scatter(-0.58, 0.34, color="grey", s=20)
    plt.show()

    x_s = np.arange(0.6, 1.6, 1.0 / k)
    y_s = np.arange(0.6, 1.6, 1.0 / k)

    x, y = np.meshgrid(x_s, y_s)
    plt.contourf(x, y, z2, levels=levels, colors=colors)
    plt.title("Root analysis")
    plt.colorbar()
    plt.scatter(1.05, 1.11, color="grey", s=20)
    plt.show()
