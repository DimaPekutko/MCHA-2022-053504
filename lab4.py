import sympy
import numpy

EPS = 0.00001

m = 0.2
a = 0.9

# Подстановка в исходные уравнения для вычисления значения f(x,y)
def val1(x, y):
    return numpy.tan(x * y + m) - x
def val2(x, y):
    return a * (x ** 2) + 2 * (y ** 2) - 1

# Выраженные из исходных уравнений x = функция(x,y), y = функция(x,y)
def eqx(x, y):
    return numpy.tan(x * y + m)
def eqy(x, y):
    return numpy.sqrt((1 - a * (x ** 2)) / 2)

# Вычисление матрицы Якоби
def W(x, y):
    return numpy.array([
        [(1 + numpy.tan(x * y + m) ** 2) * y - 1, (1 + numpy.tan(x * y + m) ** 2) * x],
        [2 * a * x, 4 * y]
    ])


def SimpleSolve(x0, y0):
    iters = 0
    (x, y) = (x0, y0)
    while True:
        iters += 1
        oldx = x
        oldy = y
        x = eqx(x, y)
        y = eqy(x, y)
        if (not (numpy.isfinite(x) and numpy.isfinite(y))):
            raise RuntimeError("Sequence {x} is divergent")
        if (max(abs(x - oldx), abs(y - oldy)) < EPS):
            return (x, y, iters)

def NewtonSolve(x0, y0):
    iters = 0
    (x, y) = (x0, y0)
    while True:
        iters += 1
        w = W(x, y)
        f = numpy.array([[val1(x, y)], [val2(x, y)]])
        deltas = numpy.linalg.solve(w, -f)
        x += deltas[0][0]
        y += deltas[1][0]
        if (not (numpy.isfinite(x) and numpy.isfinite(y))):
            raise RuntimeError("Sequence {x} is divergent")
        if (max(abs(deltas)) < EPS):
            return (x, y, iters)

def main():
    (x, y) = sympy.symbols("x y")
    eq1 = sympy.tan(x * y+a) - x
    eq2 = a*x**2 + 2 * (y ** 2) - 1

    print("Equation system: ")
    print(eq1, "= 0")
    print(eq2, "= 0")
    
    print()

    (x,y,iters) = SimpleSolve(0.5,0.5)
    print("Simple solve method:")
    print("x={:.5f}, y={:.5f} (iteration count => {})".format(x,y,iters))

    print()

    (x,y,iters) = NewtonSolve(0.5,0.5)
    print("Simple solve method:")
    print("x={:.5f}, y={:.5f} (iteration count => {})".format(x,y,iters))

main()