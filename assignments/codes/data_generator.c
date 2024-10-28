// data_generator.c
#include <math.h>

// Function to generate sqrt(x - 1) values
// x_values is an array of x values, y_values is where sqrt(x - 1) results are stored
// count is the number of values
void generate_sqrt_x_minus_1(double* x_values, double* y_values, int count) {
    for (int i = 0; i < count; i++) {
        if (x_values[i] > 1) {
            y_values[i] = sqrt(x_values[i] - 1);
        } else {
            y_values[i] = 0;  // Or some other placeholder for invalid sqrt input
        }
    }
}

