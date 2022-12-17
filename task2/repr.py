import matplotlib.pyplot as plt

if __name__ == '__main__':
    with open("dots.txt") as f:
        res_args = list(map(float, f.readline().split()))
        res_vals = list(map(float, f.readline().split()))
        inter_args = list(map(float, f.readline().split()))
        inter_vals = list(map(float, f.readline().split()))
        div_knots_args = list(map(float, f.readline().split()))
        div_knots_vals = list(map(float, f.readline().split()))
        func = list(map(float, f.readline().split()))
    y_max = max(max(res_vals), max(inter_vals))
    y_min = min(min(res_vals), min(inter_vals))
    plt.xlim(res_args[0], res_args[len(res_args) - 1])
    plt.ylim(y_min - 1, y_max + 1)
    plt.plot(res_args, res_vals, color="blue", linewidth=1, label='Approximation')
    plt.plot(res_args, func, color="orange", linewidth=0.5, label='Function')
    plt.scatter(inter_args, inter_vals, color="red", s=10, marker="o", label='Knots')
    plt.scatter(div_knots_args, div_knots_vals, color="black", s=10, marker="o", label='Div knots')
    plt.legend()
    plt.grid()
    plt.show()
