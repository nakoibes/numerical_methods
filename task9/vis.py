import matplotlib.pyplot as plt

if __name__ == '__main__':
    with open("dots.txt") as f:
        x = list(map(float, f.readline().split()))
        y1_ra = list(map(float, f.readline().split()))

    plt.xlim(x[0], x[len(x) - 1])
    plt.plot(x, y1_ra, color="orange", linewidth=0.5, label='diff')
    plt.legend()
    plt.grid()
    plt.show()
