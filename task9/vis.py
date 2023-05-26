import matplotlib.pyplot as plt

if __name__ == '__main__':
    with open("dots.txt") as f:
        x = list(map(float, f.readline().split()))
        y1_ra = list(map(float, f.readline().split()))
        y1_shoo = list(map(float, f.readline().split()))
        y1_sec = list(map(float, f.readline().split()))
        y1_new = list(map(float, f.readline().split()))

    plt.xlim(x[0], x[len(x) - 1])
    plt.plot(x, y1_ra, color="orange", linewidth=0.5, label='Diff')
    plt.plot(x, y1_shoo, color="black", linewidth=0.5, label='Shoot')
    plt.legend()
    plt.grid()
    plt.show()

    plt.xlim(x[0], x[len(x) - 1])
    plt.plot(x, y1_sec, color="green", linewidth=0.5, label='Secant')
    plt.plot(x, y1_new, color="red", linewidth=0.5, label='Newton')
    plt.legend()
    plt.grid()
    plt.show()
