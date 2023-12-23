import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    with open("dots1.txt") as f:
        m = int(f.readline())
        n = int(f.readline())
        u = [list(map(float, f.readline().split())) for i in range(n)]

    h = 1 / (m - 1)
    tau = 1 / (n - 1)

    x = [i * h for i in range(m)]
    t = [i * tau for i in range(n)]

    x, y = np.meshgrid(x, t)

    levels = [i * 0.1 for i in range(-10, 11)]

    plt.title("Explicit scheme")
    plt.contourf(x, y, u, levels=levels)
    plt.colorbar()
    plt.show()

    with open("dots2.txt") as f:
        m = int(f.readline())
        n = int(f.readline())
        u = [list(map(float, f.readline().split())) for i in range(n)]

    h = 1 / (m - 1)
    tau = 1 / (n - 1)

    x = [i * h for i in range(m)]
    t = [i * tau for i in range(n)]

    x, y = np.meshgrid(x, t)

    levels = [i * 0.1 for i in range(-10, 11)]

    plt.title("Weight scheme")
    plt.contourf(x, y, u, levels=levels)
    plt.colorbar()
    plt.show()
