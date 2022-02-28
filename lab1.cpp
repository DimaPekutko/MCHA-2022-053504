#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

const int size = 5;

double A[size][size];
double C [size][size] = {
    {0.2, 0.0, 0.2, 0.0, 0.2},
    {0.0, 0.2, 0.0, 0.2, 0.0},
    {0.2, 0.0, 0.2, .0, 0.2},
    {0.0, 0.2, 0.0, 0.2, 0.0},
    {0.0, 0.0, 0.2, 0.0, 0.2},
};
double D [size][size] = {
    {2.33, 0.81, 0.67, 0.92, -0.53},
    {-0.53, 2.33, 0.81, 0.67, 0.92},
    {0.92, -0.53, 2.33, 0.81, 0.67},
    {0.67, 0.92, -0.53, 2.33, 0.81},
    {0.81, 0.67, 0.92, -0.53, 2.33},
};
double B [size] = { 4.2, 4.2, 4.2, 4.2, 4.2 };

double X [size];
double K = 5, tmp = 0.0;

void print_A(string msg = "") {
    cout << msg << endl;
    for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			 printf("%f ", A[i][j]);
		}
        printf("| %f", B[i]);
        printf("\n");
	}
}

void gen_A() {
    for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			A[i][j] = K*C[i][j] + D[i][j];
		}
        B[i] = 4.2;
	}
}

void find_solutions() {
    for (int i = size - 1; i >= 0; i--) {
		for (int j = i; j < size-1; j++) {
			B[i] -= A[i][j+1] * X[j+1];
		}
		X[i] = B[i] / A[i][i];
	}
    printf("\nResult X vector:\n");
    for (int i = 0; i < size; i++) {
        cout << X[i] << " ";
    }
}

void just_gauss() {
    // generate A matrix
    gen_A();
    print_A("A");
    // make it triangle
    for (int i = 0; i < size; i++) {
		for (int j = i + 1; j < size; j++) {
			tmp = A[j][i] / A[i][i];
			for (int k = i + 1; k < size; k++) {
				A[j][k] -= tmp * A[i][k];
			}
			B[j] -= tmp * B[i];
			A[j][i] = 0;
		}
	}
    print_A("\nA with triangle view");
    find_solutions();
}

void gauss_partial_select() {
    gen_A();
    print_A("A");
	for (int i = 0; i < size; i++) {
		int max_index = i;
		double max = A[i][i];
		for (int j = i + 1; j < size; j++) {
			if (fabs(max) < fabs(A[j][i])) {
				max_index = j;
				max = A[j][i];
			}
		}
		// swap strings
		if (i != max_index) {
			double temp = 0;
			double temp_root = B[i];
            B[i] = B[max_index];
			B[max_index] = temp_root;
			for (int j = i; j < size; j++) {
				temp = A[i][j];
				A[i][j] = A[max_index][j];
				A[max_index][j] = temp;
			}
		}
		// to triangle vew
		for (int j = i + 1; j < size; j++) {
			tmp = A[j][i]/A[i][i];
			B[j] -= tmp * B[i];
			for (int l = i + 1; l < size; l++) {
				A[j][l] -= tmp * A[i][l];
			}
			A[j][i] = 0;
		}
	}
    print_A("\nA with triangle view");
    find_solutions();
}

void gauss_full_select() {
    gen_A();
    print_A("A");
    for (int i = 0; i < size; i++) {
		int max_index = i;
		double max = A[i][i];
		for (int j = i; j < size; j++) {
			for (int k = i; k < size; k++) {
				if (fabs(max) < fabs(A[j][k])) {
					max_index = j;
					max = A[j][k];
				}
			}
		}
		// swap strings
		if (i != max_index) {
			double temp_root = B[i];
			double temp = 0;
            B[i] = B[max_index];
			B[max_index] = temp_root;
			for (int j = i; j < size; j++) {
				temp = A[i][j];
				A[i][j] = A[max_index][j];
				A[max_index][j] = temp;
			}
		}
		// to triangle view
		for (int j = i + 1; j < size; j++) {
			tmp = A[j][i] / A[i][i];
			B[j] -= tmp * B[i];
			for (int l = i + 1; l < size; l++) {
				A[j][l] -= tmp * A[i][l];
			}
			A[j][i] = 0;
		}
	}
    print_A("\nA with triangle view");
    find_solutions();
}
 

int main()
{
	cout << fixed << setprecision(4);
    // Solving systems of linear equations by the Gauss method
    cout << "\tJust Gauss method: \n";
    	just_gauss();
    cout << "\n\n\tGauss partial select method: \n";
		gauss_partial_select();
    cout << "\n\n\tGauss full select method: \n";
		gauss_full_select();
    return 0;
}
