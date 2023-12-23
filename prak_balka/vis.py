import matplotlib.pyplot as plt

dots_sage = []
w_x_sage = []
th_x_sage = []
M_x_sage = []
Q_x_sage = []


with open("dots.txt") as f:
    Pk = int(f.readline())
    R1k = int(f.readline())
    R2k = int(f.readline())
    R3k = int(f.readline())
    q1k = int(f.readline())
    q2k = int(f.readline())
    qlenk = int(f.readline())
    Mk = int(f.readline())
    kk = int(f.readline())
    sealk = int(f.readline())
    l = float(f.readline())
    
Px = Pk*l
R1x = R1k*l
R2x = R2k*l
R3x = R3k*l
q1x = q1k*l
q2x = q2k*l
Mx = Mk*l
kx = kk*l

with open("sage_output.txt") as f:
    for line in f:
        dots_sage.append(float(line.split(',')[0]))
        w_x_sage.append(float(line.split(',')[1]))
        th_x_sage.append(float(line.split(',')[2]))
        M_x_sage.append(float(line.split(',')[3]))
        Q_x_sage.append(float(line.split(',')[4]))


dots_cpp = []
w_x_cpp = []
th_x_cpp = []
M_x_cpp = []
Q_x_cpp = []

with open("cpp_output.txt") as f:
    for line in f:
        dots_cpp.append(float(line.split(',')[0]))
        w_x_cpp.append(float(line.split(',')[1]))
        th_x_cpp.append(float(line.split(',')[2]))
        M_x_cpp.append(float(line.split(',')[3]))
        Q_x_cpp.append(float(line.split(',')[4]))

plt.subplot(5, 1, 1)
plt.xlim(0, 25*l)
plt.ylim(-1, 1)
plt.text(0, 0, 's')
plt.scatter(0, 0, color='black')
plt.text(Px, 0, 'P')
plt.scatter(Px, 0, color='green')
plt.text(R1x, 0, 'R1')
plt.scatter(R1x, 0, color='orange')
plt.text(q1x, 0, 'q1', )
plt.scatter(q1x, 0, color='yellow')
plt.text(Mx, 0, 'M')
plt.scatter(Mx, 0, color='b')
plt.text(R2x, 0, 'R2')
plt.scatter(R2x, 0, color='orange')
plt.text(q2x, 0, 'q2')
plt.scatter(q2x, 0, color='yellow')
plt.text(R3x, 0, 'R3')
plt.scatter(R3x, 0, color='orange')
plt.text(kx, 0, 'k')
plt.scatter(kx, 0, color='purple')



plt.grid(visible=True, which='major', axis='both', linewidth=0.5, color='k')
plt.grid(visible=True, which='minor', axis='both', linestyle=':')

plt.minorticks_on()
plt.subplot(5, 1, 2)
plt.xlim(0, 25*l)
plt.plot(dots_sage, Q_x_sage, color='r', linewidth=1, label='sage')
plt.plot(dots_cpp, Q_x_cpp, label='cpp')
plt.minorticks_on()
plt.ylabel('Q(x)')
plt.grid(visible=True, which='major', axis='both', linewidth=1, color='k')
plt.grid(visible=True, which='minor', axis='both', linestyle=':')

plt.subplot(5, 1, 3)
plt.xlim(0, 25*l)
plt.plot(dots_sage, M_x_sage, color='r', linewidth=1, label='sage')
plt.plot(dots_cpp, M_x_cpp, label='cpp')
plt.minorticks_on()
plt.ylabel('M(x)')
plt.grid(visible=True, which='major', axis='both', linewidth=1, color='k')
plt.grid(visible=True, which='minor', axis='both', linestyle=':')

plt.subplot(5, 1, 4)
plt.xlim(0, 25*l)
plt.plot(dots_sage, th_x_sage, color='r', linewidth=1, label='sage')
plt.plot(dots_cpp, th_x_cpp, label='cpp')
plt.ylabel('th(x)')
plt.minorticks_on()
plt.grid(visible=True, which='major', axis='both', linewidth=1, color='k')
plt.grid(visible=True, which='minor', axis='both', linestyle=':')

plt.subplot(5, 1, 5)
plt.xlim(0, 25*l)
plt.plot(dots_sage, w_x_sage, color='r', linewidth=1, label='sage')
plt.plot(dots_cpp, w_x_cpp, label='cpp')
plt.minorticks_on()
plt.ylabel('W(x)')
plt.grid(visible=True, which='major', axis='both', linewidth=1, color='k')
plt.grid(visible=True, which='minor', axis='both', linestyle=':')

for i in range(2, 6):
    plt.subplot(5, 1, i)
    plt.scatter(0, 0, color='black')
    plt.scatter(Px, 0, color='green')
    plt.scatter(R1x, 0, color='orange')
    plt.scatter(q1x, 0, color='yellow')
    plt.scatter(Mx, 0, color='b')
    plt.scatter(R2x, 0, color='orange')
    plt.scatter(q2x, 0, color='yellow')
    plt.scatter(R3x, 0, color='orange')
    plt.scatter(kx, 0, color='purple')

# plt.legend()
plt.show()
