#include <stdio.h>
#include <math.h>

#define MAX_ITER 1000
#define TOLERANCE 1e-10

void matrix_multi(int n, double A[n][n], double B[n][n], double C[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                C[i][j] += (A[i][k] * B[k][j]);
            }
        }
    }
}

double dot_product_columns(int n, double A[n][n], double Q[n][n], int i, int j) {
    double result = 0.0;
    for (int k = 0; k < n; k++) {
        result += A[k][i] * Q[k][j];
    }
    return result;
}

void qr_decomposition(int n, double A[n][n], double Q[n][n], double R[n][n]) {
    double u[n][n];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Q[i][j] = 0.0;
            R[i][j] = 0.0;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            u[j][i] = A[j][i];
        }

        for (int k = 0; k < i; k++) {
            R[k][i] = dot_product_columns(n, A, Q, i, k);
            for (int j = 0; j < n; j++) {
                u[j][i] -= R[k][i] * Q[j][k];
            }
        }
        double norm = sqrt(dot_product_columns(n, u, u, i, i));
	if(norm<TOLERANCE){
		return;
	}
        R[i][i] = norm;

        for (int j = 0; j < n; j++) {
            Q[j][i] = u[j][i] / norm;
        }
    }
}

void find_eigenvalues(int n, double A[n][n]) {
    double Q[n][n], R[n][n], A1[n][n];
    int iter = 0;

    while (iter < MAX_ITER) {
        double shift = A[n - 1][n - 1];
        for (int i = 0; i < n; i++) {
            A[i][i] -= shift;
        }

        qr_decomposition(n, A, Q, R);

        matrix_multi(n, R, Q, A1);

        double off_diagonal_sum = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    off_diagonal_sum += fabs(A[i][j]);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = (i == j ? A1[i][j] + shift : A1[i][j]);
            }
	}

        if (off_diagonal_sum < TOLERANCE) {
            break;
        }
        iter++;
    }

    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("%.10lf\n", A[i][i]);
    }
}

int main() {
    int n;
    printf("Enter the order of the matrix: ");
    scanf("%d", &n);

    if (n <= 0) {
        printf("Invalid matrix size.\n");
        return 1;
    }

    double A[n][n];
    printf("Enter the elements of the matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &A[i][j]);
        }
    }

    find_eigenvalues(n, A);
    return 0;
}

