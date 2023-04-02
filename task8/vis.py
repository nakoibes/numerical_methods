import matplotlib.pyplot as plt

if __name__ == '__main__':
    with open("dots.txt") as f:
        f_args = list(map(float, f.readline().split()))
        f_vals = list(map(float, f.readline().split()))


    # y_max = max(max(f_vals), max(a1_vals),max(a2_vals),max(ar_vals))
    # y_min = min(min(f_vals), min(a1_vals),min(a2_vals),min(ar_vals))
    plt.xlim(f_args[0], f_args[len(f_args) - 1])
    # plt.ylim(y_min - 1, y_max + 1)
    # plt.plot(f_args, f_vals, color="blue", linewidth=1, label='dy/dx')
    plt.plot(f_args, f_vals, color="orange", linewidth=0.5, label='Function')
    # plt.scatter(a1_args, a1_vals, color="red", s=10, marker="o", label='Step h')
    # plt.scatter(a2_args, a2_vals, color="green", s=10, marker="o", label='Step h/2')
    # plt.scatter(ar_args, ar_vals, color="black", s=10, marker="o", label='Runge')
    # plt.scatter(div_knots, [0]*len(div_knots), color="green", s=15, marker="o", label='Div_knots')
    plt.legend()
    plt.grid()
    plt.show()