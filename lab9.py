import math
from scipy.integrate import odeint
import numpy
from matplotlib import pylab
a = 0
b = 1
y0 = 0
epsilon = 0.0001


def f(y, x):
    return (0.9*(1-y**2))/((1+1.5)*(x**2)+(y**2)+1)

def runge_kutta(x, step):
    if x <= a:
        return y0
    prev_x = x - step
    prev_y = runge_kutta(prev_x, step)
    F1 = f(prev_y, prev_x)
    F2 = f(prev_y + (step / 2) * F1, prev_x + step / 2)
    F3 = f(prev_y + (step / 2) * F2, prev_x + step / 2)
    F4 = f(prev_y + step * F3, x)
    return prev_y + step * (F1 + 2 * F2 + 2 * F3 + F4) / 6
    
def find_step(epsilon):
    h0 = epsilon ** (1/4)
    n = int((b - a) // h0)
    if n % 2 != 0 :
        n += 1
    while check_step(n, epsilon):
        n = n // 4 * 2
    while not check_step(n, epsilon):
        n += 2
    return (b - a) / n

def check_step(n, epsilon):
    h = (b - a) / n     
    y2 = runge_kutta(a + 2 * h, h)
    y2e = runge_kutta(a + 2 * h, h * 2)
    eps = (1/15) * abs(y2 - y2e)
    return eps < epsilon

step = find_step(epsilon)
print("РЁР°Рі РёРЅС‚РµРіСЂРёСЂРѕРІР°РЅРёСЏ: ", step)

def accurate_solution(x):
    sol = odeint(f, y0, [a, x])
    return sol[1][0]

def show_curve(method, step, text):
    pylab.cla()
    xlist = numpy.arange(a, b + step, step)
    ylist = []
    x = a
    while x <= b:
        r = method(x, step)
        ylist.append(r)
        x += step
    pylab.plot (xlist, ylist, label = text)
    print(f"РњРµС‚РѕРґРѕРј Р СѓРЅРіРµ-РљСѓС‚С‚Р°: {ylist}")
    pylab.grid(True)
    pylab.legend()
    pylab.savefig("1.png")
    ylist2 = ylist
    ylist = []
    x = a
    while x <= b:
        r = accurate_solution(x)
        ylist.append(r)
        x += step
    pylab.cla()
    pylab.plot (xlist, ylist, label = "РўРѕС‡РЅРѕРµ СЂРµС€РµРЅРёРµ")
    print(f"РўРѕС‡РЅРѕРµ СЂРµС€РµРЅРёРµ: {ylist}")
    pylab.grid(True)
    pylab.legend()
    pylab.savefig("2.png")
    rez = []
"""    for item in range(len(ylist)):
        rez.append(abs(ylist[item]-ylist2[item]))
    print(f"РџРѕРіСЂРµС€РЅРѕСЃС‚СЊ РјРµР¶РґСѓ Р·Р°С‡РµРЅРёСЏРјРё: {rez}")"""

show_curve(runge_kutta, step, "РєСЂРёРІР°СЏ РјРµС‚РѕРґРѕРј Р СѓРЅРіРµ-РљСѓС‚С‚Р°")

def euler(x, step):
    if x <= a:
        return y0
    prev = euler(x - step, step)
    return prev + step * f(prev, x)

xlist = numpy.arange(a, b + step, step)
runge_kutta_points = []
euler_points = []
accurate_solution_points = []
x = a
while x <= b:
    r = runge_kutta(x, step)
    r1 = euler(x, step)
    r2 = accurate_solution(x)
    runge_kutta_points.append(r)
    euler_points.append(r1)
    accurate_solution_points.append(r2)
    x += step

print(f"Euler points: {euler_points}")
"""rez = []
for item in range(len(euler_points)):
    rez.append(abs(euler_points[item]-runge_kutta_points[item]))"""

print(f"Result: {rez}")
pylab.cla()
pylab.plot (xlist, runge_kutta_points,label = "", color = (0, 0, 1))
pylab.plot (xlist, euler_points,label = "", color = (1, 0, 0))
pylab.plot (xlist, accurate_solution_points,label = "", color = (0, 1, 0))
pylab.grid(True)
pylab.legend()
pylab.savefig("combine.png")
