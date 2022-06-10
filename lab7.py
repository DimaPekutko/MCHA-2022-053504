# import numpy as np
# import math
# import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import math

dataset = []
def spline_interpolation(data, _x_, tol = 1e-100):
    if not False:
        x = np.array(data[0])
        y = np.array(data[1])
        if np.any(np.diff(x) < 0):
            idx = np.argsort(x)
            x = x[idx]
            y = y[idx]
        _range_ = len(x)
        dx = np.diff(x)
        dy = np.diff(y)
        A = np.zeros(shape=(_range_,_range_))
        b = np.zeros(shape=(_range_,1))
        A[0,0] = 1
        A[_range_-1][_range_-1] = 1
        for i in range(1,_range_-1):
            A[i, i-1] = dx[i-1]
            A[i, i+1] = dx[i]
            A[i,i] = 2*(dx[i-1]+dx[i])
            b[i,0] = 3*(dy[i]/dx[i] - dy[i-1]/dx[i-1])
        c = np.linalg.solve(A, b)
        d = np.zeros(shape = (_range_-1,1))
        b = np.zeros(shape = (_range_-1,1))
        for i in range(0,len(d)):
            d[i] = (c[i+1] - c[i]) / (3*dx[i])
            b[i] = (dy[i]/dx[i]) - (dx[i]/3)*(2*c[i] + c[i+1])    
        func_idx = 0
        for i in range(_range_ - 1):
            if _x_ > x[i]:
                func_idx = i
                break
    return {"result":y[i] + b[i]*(_x_-x[i]) + c[i]*(_x_-x[i])**2 + d[i]*(_x_-x[i])**3,  "coef":(b.squeeze(), c.squeeze(), d.squeeze())}
# plt.plot(x,y,'ro')
# plt.plot(x,y, 'b')
# plt.title("Data set and linear interpolation")
# plt.show()

# x = [0, 1, 2]
# y = [1, 3, 2]
def build_spline(title, x, y, a, b):
    f = interpolate.CubicSpline(x, y)
    x_new = np.linspace(a, b, 100)
    y_new = f(x_new)
    plt.figure(figsize = (5,5))
    plt.plot(x_new, y_new, 'b')
    plt.plot(x, y, 'ro')
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    print(title, "y((b-a)/2) =", round(float(f((b-a)/2)),4))
    print("x = ", x)
    print("y = ", y)
    print("__________")

dataset1 = [
    [i*(2/5) for i in range(6)],
    [math.sinh(i*(2/5)) for i in range(6)]
]

dataset2 = [
    [i for i in range(5)],
    [math.sqrt(i) for i in range(5)]
]


build_spline("Cubic spline interp for sinh(x)", dataset1[0], dataset1[1], 0, 2)
build_spline("Cubic spline interp for sqrt(x)", dataset2[0], dataset2[1], 0, 4)


# s = interpolate.InterpolatedUnivariateSpline(x, y)
# xfit = np.arange(0, n-1, np.pi/50)
# yfit = s(xfit)

# print(2)
# plt.plot(x, y, 'ro')
# plt.plot(xfit, yfit,'green')
# plt.title("InterpolatedUnivariateSpline interpolation")
# plt.show() 