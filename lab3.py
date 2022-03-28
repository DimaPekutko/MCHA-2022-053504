import sympy as sp
from sympy.abc import x


def section(func, border):
    counter = 0
    size = len(func)
    flag = func[0](border) > 0
    for i in range(size - 1):
        change_flag = func[i + 1](border) > 0
        if change_flag != flag:
            counter += 1
        flag = change_flag
    return counter


def shturm(polynomial, left_border, right_border):
    func_x = []
    func_x.append(polynomial)
    func_x.append(sp.diff(polynomial))
    degree = sp.degree(polynomial, gen=x)
    for i in range(degree - 1):
        div = sp.div(func_x[i], func_x[i + 1])
        func_x.append(-div[1])
    amount = section(func_x, left_border) - section(func_x, right_border)
    return amount


def bisection(polynomial, left_border, right_border, eps):
    iteration = 0
    diff = 10 * eps
    mid = 0
    base_iter = []
    base_diff = []
    while diff > eps:
        mid = (left_border + right_border) / 2
        if polynomial(left_border) * polynomial(mid) <= 0:
            right_border = mid
        else:
            left_border = mid
        diff = abs((left_border - right_border))
        base_iter.append(iteration)
        base_diff.append(diff)
        iteration += 1
    # plot_graph_scatter(base_iter, base_diff)
    return mid, iteration


def chord(polynomial, left_border, right_border, eps):
    iteration = 0
    diff = 10 * eps
    mid = 10 * eps
    base_iter = []
    base_diff = []
    while diff > eps:
        temp = mid
        mid = left_border - (polynomial(left_border) / (polynomial(right_border) - polynomial(left_border))) * (
            right_border - left_border)
        if polynomial(left_border) * polynomial(mid) <= 0:
            right_border = mid
        else:
            left_border = mid
        diff = abs((mid - temp))
        base_iter.append(iteration)
        base_diff.append(diff)
        iteration += 1
    return mid, iteration


def newton(polynomial, border, eps):
    iteration = 0
    diff = 10 * eps
    base_iter = []
    base_diff = []
    while diff > eps:
        border_temp = border
        border = border - polynomial(border) / sp.diff(polynomial)(border)
        diff = abs(border - border_temp)
        base_iter.append(iteration)
        base_diff.append(diff)
        iteration += 1
    return border, iteration


def main():
    a = sp.Float(9.57496)
    b = sp.Float(-243.672)
    c = sp.Float(773.65)
    border = [-10, 10]
    eps = 0.0001
    polynomial = sp.poly(x ** 3 + a * x ** 2 + b * x + c)
    res = sp.solve_poly_system([polynomial])
    print(res)
    print("Roots count: ", shturm(polynomial, border[0], border[1]))
    bis_res = bisection(polynomial, border[0], border[0] / 2, eps)
    print("Bisection: ", round(bis_res[0], 4), ": Iteration =", bis_res[1])
    chords_res = chord(polynomial, border[0], border[0] / 2, eps)
    print("Chords: ", round(chords_res[0], 4), ": Iteration =", chords_res[1])
    newton_res = newton(polynomial, border[0], eps)
    print("Bisection: ", round(newton_res[0], 4), ": Iteration =", newton_res[1])


if __name__ == "__main__":
    main()
