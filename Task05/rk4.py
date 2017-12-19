from math import sqrt

##########################################
# classical Runge–Kutta method
##########################################
def rk4(f, x0, y0, x1, n):
    """
    f - function;
    x0, y0 - initial values;
    x1 - end value;
    n - integration step;
    """
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (x1 - x0) / float(n)
    vx[0] = x = x0
    vy[0] = y = y0
    for i in range(1, n + 1):
        k1 = h * f(x, y)
        k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
        k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
        k4 = h * f(x + h, y + k3)
        vx[i] = x = x0 + i * h
        vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6
    return vx, vy

##########################################
# input function
##########################################
def f(x, y):
    return x * sqrt(y)

##########################################
# test data
##########################################
vx, vy = rk4(f, 0, 1, 10, 100)

def exact_solution(x):
    return (4 + x * x)**2 / 16

##########################################
# output
##########################################
print("x; y; error")
for x, y in list(zip(vx, vy))[::10]:
    print("{0:.2f};{1:10.2f};{2:10f}".format(x, y, y - exact_solution(x)))