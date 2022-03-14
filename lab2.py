# Method of simple iterations && Seidel's method

from math import fabs
from statistics import mean

C = [
    [0.01,0,-0.02,0,0],
    [0.01,0.01,-0.02,0,0],
    [0,0.01,0.01,0,-0.02],
    [0,0,0.01,0.01,0],
    [0,0,0,0.01,0.01]
]

D = [
    [1.33,0.21,0.17,0.12,-0.13],
    [-0.13,-1.33,0.11,0.17,0.12],
    [0.12,-0.13,-1.33,0.11,0.17],
    [0.17,0.12,-0.13,-1.33,0.11],
    [0.11,0.67,0.12,-0.13,-1.33]
]


b = [1.201,2.2,4.0,0.0,-1.2]

VARIANT = 5
SIZE = len(C[0])
PREC = 5
MAX_ITERATION_COUNT = 100000

def print_matrix(matrix):
    for i in range(SIZE):
        print(matrix[i])

def simple_iteration_method(matrix):
    x_vector = [0,0,0,0,0]
    for k in range(MAX_ITERATION_COUNT):
        buf_x_vector = [*x_vector]
        # print(buf_x_vector)
        is_acc_x = True
        for i in range(SIZE):
            str_sum = 0
            for j in range(SIZE):
                str_sum += matrix[i][j]*x_vector[j]
            buf_x_vector[i] = (b[i]-str_sum)/matrix[i][i]+x_vector[i]
            if fabs(round(buf_x_vector[i], PREC) - round(x_vector[i], PREC)):
                is_acc_x = False
        if is_acc_x:
            return [round(x, PREC) for x in x_vector], k
        x_vector = buf_x_vector
    return [0]*SIZE, k

def seidel_method(matrix):
    x_vector = [0,0,0,0,0]
    for k in range(MAX_ITERATION_COUNT):
        buf_x_vector = [*x_vector]
        is_acc_x = True
        for i in range(SIZE):
            str_sum = 0
            for j in range(SIZE):
                str_sum += matrix[i][j]*buf_x_vector[j]
            buf_x_vector[i] = (b[i]-str_sum)/matrix[i][i]+x_vector[i]
            if fabs(round(buf_x_vector[i], PREC) - round(x_vector[i], PREC)):
                is_acc_x = False
        if is_acc_x:
            return [round(x, PREC) for x in x_vector], k
        x_vector = buf_x_vector
    return [0]*SIZE, k

def main():
    standard_matrix = [
        [100,2,100,4,5],
        [10,99,8,7,6],
        [12,23,100,65,77],
        [2,1,11,96,33],
        [96,32,4,5,96]
    ]

    # # generate standard matrix
    for i in range(SIZE):
        for j in range(SIZE):
            standard_matrix[i][j] = (round(C[i][j]*VARIANT+D[i][j],PREC))
            # standard_matrix[i][j] = i+j+1
    
    # standard_matrix[3][3] = 0

    for i in range(SIZE):
        if standard_matrix[i][i] == 0.0:
            exit("There a zero at matrix diagonal.")

    # # matrix norm
    # str_sum = 0
    for i in range(SIZE):
        buffer_sum = 0
        for j in range(SIZE):
            if i != j:
                buffer_sum += fabs(standard_matrix[i][j])
        print(buffer_sum)
    
    iter_x_vector, iter_count = simple_iteration_method(standard_matrix)
    seidel_x_vector, seidel_count = seidel_method(standard_matrix)

    

    print("A matrix:")
    print_matrix(standard_matrix)
    if mean(iter_x_vector):
        print(f"Iteration method. X result vector (steps count = {iter_count}):")
        print_matrix(iter_x_vector)
    else:
        print(f"Iteration method. Cannot find solution {iter_x_vector}")
    if mean(seidel_x_vector):
        print(f"Seidel's method. X result vector (steps count = {seidel_count}):")
        print_matrix(seidel_x_vector)
    else:
        print(f"Seidel's method. Cannot find solution {iter_x_vector}")


if __name__ == "__main__":
    main()