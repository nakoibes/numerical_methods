import matplotlib.pyplot as plt

if __name__ == '__main__':
    with open("dots.txt") as f:
        x = list(map(float, f.readline().split()))
        y1_e = list(map(float, f.readline().split()))
        y1_e_m = list(map(float, f.readline().split()))
        y1_rk2 = list(map(float, f.readline().split()))
        y1_rk4 = list(map(float, f.readline().split()))
        y1_a = list(map(float, f.readline().split()))
        y2_e = list(map(float, f.readline().split()))
        y2_e_m = list(map(float, f.readline().split()))
        y2_rk2 = list(map(float, f.readline().split()))
        y2_rk4 = list(map(float, f.readline().split()))
        y2_a = list(map(float, f.readline().split()))

    plt.xlim(x[0], x[len(x) - 1])
    plt.plot(x, y1_e, color="orange", linewidth=0.5, label='eil')
    plt.plot(x, y1_e_m, color="red", linewidth=0.5, label='eil mod')
    plt.plot(x, y1_rk2, color="black", linewidth=0.5, label='rk2')
    plt.plot(x, y1_rk4, color="blue", linewidth=0.5, label='rk4')
    plt.plot(x, y1_a, color="green", linewidth=0.5, label='adamas')

    plt.legend()
    plt.grid()
    plt.show()
    plt.plot(x, y2_e, color="orange", linewidth=0.5, label='eil')
    plt.plot(x, y2_e_m, color="red", linewidth=0.5, label='eil mod')
    plt.plot(x, y2_rk2, color="black", linewidth=0.5, label='rk2')
    plt.plot(x, y2_rk4, color="blue", linewidth=0.5, label='rk4')
    plt.plot(x, y2_a, color="green", linewidth=0.5, label='adamas')
    plt.legend()
    plt.grid()
    plt.show()
